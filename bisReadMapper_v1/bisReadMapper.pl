#!/usr/bin/perl
use strict;
# bisReadMapper.pl: a perl script to map single-end bisulfite sequencing reads and report the methylation levels and SNPs. 
# USAGE: ./bisReadMapper.pl params.txt > SAMPLE_NAME.log
# This version performs trimming and the old samtools pileup command.

my $ref_dir = 0;
my @reads = ();
my @refList = ();
my $snp_file = 0;
my $align_mode = 0;
my $qual_base = 0;
my $cpu = 1;
my $read_len = 0;
my $soap_dir = 0;
my $samtools_dir = 0;
my $name = "Sample";
my $allC = 0;
my $threep = 0;
my $fivep= 0;
my $bam=0;

my $script_dir = `readlink -f $0`;
chop($script_dir);
$script_dir =~ s/\/bisReadMapper.pl//g;

my $keep_bam = 1;
my $minDepth = 10;

my ($template_fwd, $template_fwd_fa, $template_rev, $template_rev_fa, $template_idx);

my ($soap2_exe, $samtools, $soap2sam, $bcftools);

my %rcTable;
$rcTable{'A'}='T';
$rcTable{'T'}='A';
$rcTable{'G'}='C';
$rcTable{'C'}='G';
$rcTable{'N'}='N';
$rcTable{'R'}='Y';
$rcTable{'Y'}='R';
$rcTable{'M'}='K';
$rcTable{'K'}='M';
$rcTable{'S'}='S';
$rcTable{'W'}='W';

my %chrSizes;

sub main(){
	open(PARAMS, "$ARGV[0]") || die("No params file given or params file is unreadable\n");
	while(my $a = <PARAMS>){
		next if($a =~ m/#/);
		chop($a);
		if($a =~ m/reads=/){ $a =~ s/reads=//g; @reads = split(",", $a); print "Reads: ", join(",", @reads), "\n";}
		if($a =~ m/length=/){ $a =~ s/length=//g; $read_len = int($a); print "Read lengths: $read_len\n"; }
		if($a =~ m/refDir=/){ $a =~ s/refDir=//g; $ref_dir = $a; print "Reference dir: $ref_dir\n";}
		if($a =~ m/soapDir=/){ $a =~ s/soapDir=//g; $soap_dir = $a; print "SOAP dir: $soap_dir\n";}
		if($a =~ m/alignMode=/){ $a =~ s/alignMode=//g; $align_mode = $a; print "Reads alignment mode: $align_mode\n";}
		if($a =~ m/qualBase=/) { $a =~ s/qualBase=//g; $qual_base = $a; print "Base quality: $qual_base\n";}
		if($a =~ m/numCPU=/){ $a =~ s/numCPU=//g; $cpu = int($a); print "Number of processors for mapping: $cpu\n";}
		if($a =~ m/samtoolsDir=/){ $a =~ s/samtoolsDir=//g; $samtools_dir = $a; print "Samtools dir: $samtools_dir\n";}
		if($a =~ m/soap2sam=/)   { $a =~ s/soap2sam=//g; $soap2sam = $a; print "Soap2sam path: $soap2sam\n";}
		if($a =~ m/snp=/){ $a =~ s/snp=//g; $snp_file = $a; print "SNP file= $snp_file\n";}
		if($a =~ m/name=/){ $a =~ s/name=//g; $name = $a; print "Sample name= $name\n";}
		if($a =~ m/allC=/){ $a =~ s/allC=//g; $allC = int($a); print "Call all C ? = $allC\n";}
		if($a =~ m/trim3=/){ $a =~ s/trim3=//g; $threep = int($a); print "Trim 3' = $threep\n";}
		if($a =~ m/trim5=/){ $a =~ s/trim5=//g; $fivep = int($a); print "Trim 5' = $fivep\n";}
		if($a =~ m/bam=/){ $a =~ s/bam=//g; $bam = $a; print "Input is BAM? = $bam\n";}
	}
	close(PARAMS);
	if(!$align_mode){
		die("Need to specify alignment type.\n");
	}
	if(!$qual_base){
		die("Need to specify base quality scale!\n");
	}
	if(!$reads[0]){
		die("No reads provided.\n");
	}	
	if($align_mode eq 'P' && scalar(@reads) < 2){
		die("Can't align by paired-end, only one reads file provided.\n");
	}
	if(!$ref_dir){
		die("Missing path to reference.\n");
	}
	if(!$soap_dir){
		die("Missing path to SOAP.\n");
	}
	if(!$samtools_dir){
		die("Missing path to samtools.\n");
	}

	$soap2_exe = $soap_dir. "/soap";
	$samtools = $samtools_dir . "/samtools";

	my $cmd = "ls $ref_dir";
	my $files = `$cmd`;
	@refList = split(/[\t\n]/, $files);

	foreach my $f (@refList){
		if($f =~ m/.bis.fwd.index.ann/g){
			$f=~s/.index.ann//g;
			$template_fwd_fa = $ref_dir."/".$f;
			print $template_fwd_fa,"\n";
		}
		if($f =~ m/.bis.rev.index.ann/g){
			$f =~ s/.index.ann//g;
			$template_rev_fa = $ref_dir."/".$f;
			print $template_rev_fa,"\n";
		}
	}

	die("Missing reference files! \n Did you run genomePrep.pl?\n Check that given reference directory is correct!\n") 
	if(!$template_fwd_fa || !$template_rev_fa );

	$template_fwd = $template_fwd_fa.".index";
	$template_rev = $template_rev_fa.".index";
	$template_idx = $template_fwd_fa.".fai";

	open(GENOME_INDEX, "$template_idx") || die("Error opening chromosome sizes file for reading!\n");
	while(my $line = <GENOME_INDEX>){
		chop($line);
		my @f = split "\t", $line;
		$chrSizes{$f[0]} = $f[1];
	}
	close(GENOME_INDEX);

	if($bam == 0){
		my ($soap_fwd_map_file, $soap_rev_map_file) = (0,0);
		if($align_mode eq "P"){
			($soap_fwd_map_file, $soap_rev_map_file) = fastq2SOAPpe($reads[0], $reads[1]);
		}else{
			($soap_fwd_map_file, $soap_rev_map_file) = fastq2SOAPse($reads[0]);
		}
		extractCG_SNP($soap_fwd_map_file, $soap_rev_map_file);
	}else{
		extractCG_SNP($reads[0], $reads[1]);
	}
	if($snp_file){
		#filter SNPs
		print "Filtering SNPs\n";
		$cmd = "$script_dir/bisSnpFilter.pl $name.snp $snp_file > $name.snp.filtered";
		print $cmd, "\n";
		system($cmd) == 0 or die "system problem (exit $?): $!\n";

	}
	
	#methylFreq2BED	
	print "Creating BED files for CG, CHG, CHH positions separately\n";
	$cmd = "less $name.*.methylFreq | $script_dir/methylFreq2BED.pl $name $minDepth";
	print $cmd, "\n";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
		
	print "Calculating bisulfite data quality by correlating the forward and reverse methylation frequencies for CpGs at $minDepth minimum depth \n";
	$cmd = "less $name.*.methylFreq | $script_dir/frMethylCorr.pl $minDepth";
	print $cmd, "\n";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
	
}

sub extractCG_SNP(){
	my ($fwd_map, $rev_map) = @_;
	if($bam == 1){
		$sorted_fwd_bam = $fwd_map;
		$sorted_rev_bam = $rev_map;
		my $cmd = "$samtools index $sorted_fwd_bam";
		system($cmd) == 0 or die "system problem (exit $?): $!\n";
		$cmd = "$samtools index $sorted_rev_bam";
		system($cmd) == 0 or die "system problem (exit $?): $!\n";
		my $vcf_fwd = $sorted_fwd_bam . ".vcf";
		my $vcf_rev = $sorted_rev_bam . ".vcf";
		$cmd = "$samtools pileup -cf $template_fwd_fa $sorted_fwd_bam > $vcf_fwd";
		print $cmd, "\n";
		system($cmd) == 0 or die "system problem (exit $?): $!\n";
		$cmd = "$samtools pileup -cf $template_rev_fa $sorted_rev_bam > $vcf_rev";
		print $cmd, "\n";
		system($cmd) == 0 or die "system problem (exit $?): $!\n";

		my $snpcall_file = $name . ".snp";
		open( my $snp_h, ">$snpcall_file") || die("Error writing to snp file, $snp_file \n");	
		foreach my $val (@refList){
			if($val ~ /positions.txt/){
				my $methylFreq_file = $name . "." . $val . ".methylFreq";	
				open( my $cpg_h, ">$methylFreq_file") || die("Error writing to methylation frequency file, $methylFreq_file \n");
				extractbychr($vcf_fwd, $vcf_rev, $val, $snp_h, $cpg_h);
				close($cpg_h);
			}
		}
		extractbyChr($vcf_fwd, $vcf_rev, $val);
		close($snp_h);
	
		unlink($sorted_fwd_bam . ".bai");
		unlink($sorted_rev_bam . ".bai");

		unlink($sorted_fwd_bam . ".vcf");
		unlink($sorted_rev_bam . ".vcf");
		return;
	}
	my $fqName = $fwd_map;
	$fqName =~ s/.soap.out//g;
	#merge the two soap mapped files
	open(SOAP_REV, "$rev_map") || die("Error in opening $rev_map.");
	open(SOAP_FWD, ">>$fwd_map") || die("Error in opening $fwd_map.");
	while(my $line = <SOAP_REV>){
		chop($line);
		print SOAP_FWD $line, "\tREV\n";
	}
	close(SOAP_REV);
	close(SOAP_FWD);	

	#identify unique reads and save them in two files based on the templates.
	my $soap_combined_sorted_map_file = $fqName.".combined.sorted.soap.out";
	my $cmd = "sort -k 1,1 < $fwd_map > $soap_combined_sorted_map_file";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";

	my $soap_clean_fwd_map_file = $fqName.".clean.soap.fwd.out";
	my $soap_clean_rev_map_file = $fqName.".clean.soap.rev.out";
	open(SOAP_OUT, "$soap_combined_sorted_map_file") || die("Error in opening $soap_combined_sorted_map_file.");
	open(CLEAN_SOAP_FWD, ">$soap_clean_fwd_map_file") || die("Error in opening $soap_clean_fwd_map_file.");
	open(CLEAN_SOAP_REV, ">$soap_clean_rev_map_file") || die("Error in opening $soap_clean_rev_map_file.");
	my $last_line = <SOAP_OUT>;
	my @last_fields = split(/\t/, $last_line);
	while(my $line =  <SOAP_OUT>){
		my @fields = split(/\t/, $line);
		if(!$last_line){
			$last_line = $line;
			@last_fields = @fields;
			next;
		}
		if($fields[0] eq $last_fields[0]){
			if(substr($fields[9],0,1) < substr($last_fields[9],0,1)){
				#2nd line is a better hit
				$last_line = $line;
				@last_fields = @fields;
			}elsif(substr($fields[9],0,1) == substr($last_fields[9],0,1)){
				#two equivalent hits, remove the reads
				undef($last_line);
				undef(@last_fields);
			}else{
				#1st line is a better hit, do nothing
			}
		}else{			
			my ($id,$orig_seq) = split(/\|/, $last_fields[0]);
			$last_fields[0] = $id;
			$last_fields[1] = $orig_seq;
			my $tag = pop(@last_fields);
			if($tag eq "REV\n"){
				print CLEAN_SOAP_REV join("\t", @last_fields), "\n";
			}else{
				print CLEAN_SOAP_FWD join("\t", @last_fields), "\t", $tag;
			}
			$last_line = $line;
			@last_fields = @fields;
		}
	}
	#print the last line
	if($last_line){
		my ($id, $orig_seq) = split (/\|/, $last_fields[0]);
		$last_fields[0] = $id;
		$last_fields[1] = $orig_seq;
		my $tag = pop(@last_fields);
		if($tag eq "REV\n"){
			print CLEAN_SOAP_REV join("\t", @last_fields), "\n";
		}else{
			print CLEAN_SOAP_FWD join("\t", @last_fields), "\t", $tag;
		}
	}
	close(SOAP_OUT);
	close(CLEAN_SOAP_FWD);
	close(CLEAN_SOAP_REV);
	unlink($soap_combined_sorted_map_file);
	unlink($fwd_map);
	unlink($rev_map);
	
	my $sorted_fwd_bam = $fwd_map . ".sorted";
	my $sorted_rev_bam = $rev_map . ".sorted";
	
	print $template_idx, "\n";

	my $cmd;
	if($align_mode eq 'P'){
		$cmd = "$soap2sam p $soap_clean_fwd_map_file | $samtools view -uSt $template_idx - | $samtools sort - $sorted_fwd_bam";
		print $cmd, "\n";
		system($cmd) == 0 or die "system problem (exit $?): $!\n";

		$cmd = "$soap2sam p $soap_clean_rev_map_file | $samtools view -uSt $template_idx - | $samtools sort - $sorted_rev_bam";
		print $cmd, "\n";
		system($cmd) == 0 or die "system problem (exit $?): $!\n";
	}else{
		$cmd = "$soap2sam $soap_clean_fwd_map_file | $samtools view -uSt $template_idx - | $samtools sort - $sorted_fwd_bam";
		print $cmd, "\n";
		system($cmd) == 0 or die "system problem (exit $?): $!\n";

		$cmd = "$soap2sam $soap_clean_rev_map_file | $samtools view -uSt $template_idx - | $samtools sort - $sorted_rev_bam";
		print $cmd, "\n";
		system($cmd) == 0 or die "system problem (exit $?): $!\n";
	}

	#die("Test samtools\n");
	unlink($soap_clean_fwd_map_file);
	unlink($soap_clean_rev_map_file);

	$sorted_fwd_bam = $sorted_fwd_bam . ".bam";
	$sorted_rev_bam = $sorted_rev_bam . ".bam";

	$cmd = "$samtools index $sorted_fwd_bam";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
	$cmd = "$samtools index $sorted_rev_bam";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
	my $vcf_fwd = $sorted_fwd_bam . ".vcf";
	my $vcf_rev = $sorted_rev_bam . ".vcf";
	$cmd = "$samtools pileup -cf $template_fwd_fa $sorted_fwd_bam > $vcf_fwd";
	print $cmd, "\n";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
	$cmd = "$samtools pileup -cf $template_rev_fa $sorted_rev_bam > $vcf_rev";
	print $cmd, "\n";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";

	my $snpcall_file = $name . ".snp";
	open( my $snp_h, ">$snpcall_file") || die("Error writing to snp file, $snp_file \n");	
	foreach my $val (@refList){
			if($val ~ /positions.txt/){
				my $methylFreq_file = $name . "." . $val . ".methylFreq";	
				open( my $cpg_h, ">$methylFreq_file") || die("Error writing to methylation frequency file, $methylFreq_file \n");
				extractbychr($vcf_fwd, $vcf_rev, $val, $snp_h, $cpg_h);
				close($cpg_h);
			}
	}
	extractbyChr($vcf_fwd, $vcf_rev, $val);
	close($snp_h);
	
	unlink($sorted_fwd_bam . ".bai");
	unlink($sorted_rev_bam . ".bai");

	unlink($sorted_fwd_bam . ".vcf");
	unlink($sorted_rev_bam . ".vcf");
}

sub extractbyChr(){
	my ($vcf_fwd, $vcf_rev, $val, $snp_h, $cpg_h) = @_;
	my %cpgTable;
	my $chr = "NA";
	my $pos = 0;
	my $fwd_line;
	my $rev_line;
	my $cpg_line;
	open(FWD_FILE, "$vcf_fwd") || die("Error opening file $vcf_fwd.\n");
	open(REV_FILE, "$vcf_rev") || die("Error opening file $vcf_rev.\n");
	open(C_FILE, "$ref_dir/$val") || die("Error opening filr $ref_dir/$val.\n");
	while($cpg_line = <C_FILE>){
		chop($cpg_line);
		my @fields = split(/\t/, $cpg_line);
		$pos = $fields[1];
		my $context = "CG";
		$context = $fields[2] if($fields[2]);
		next if($context ne "CG" && $allC == 0);
		my @tmp = split(/:/, $fields[0]);
		$chr = $tmp[0];
		$cpgTable{$fields[0]."\t".$pos}=$context;
	}
	#extract CpG/SNP from the fwd template	
	my $switch = 0;
	while($fwd_line = <FWD_FILE>){
		chop($fwd_line);
		my @fields = split /\t/, $fwd_line;
		$switch=1 if($fields[0] eq $chr);
		last if($fields[0] ne $chr && $switch == 1); #already finished the current chromosome;
		next if($fields[0] ne $chr);
		makeCalls("W", \%cpgTable, $fwd_line, $cpg_h, $snp_h);
	}
	#extract CpG/SNP from the rev template
	$switch = 0;
	while($rev_line = <REV_FILE>){
		chop($rev_line);
		my @fields = split /\t/, $rev_line;
		$switch=1 if($fields[0] eq $chr);
		last if($fields[0] ne $chr && $switch == 1); #already finished the current chromosome;
		next if($fields[0] ne $chr);
		makeCalls("C", \%cpgTable, $rev_line, $cpg_h, $snp_h);
	}
	undef %cpgTable;
	close(REV_FILE);
	close(FWD_FILE);
}
				
sub fastq2SOAPse(){
	my $tmp = shift;
	my @f = split /\//, $tmp;		
	my $fqName = pop(@f);
	my $encodedFqName = $fqName . ".encoded";

	encodeFastq($tmp, $encodedFqName);

	my $maxMismatches=1;
	$maxMismatches = int($read_len/40) if(int($read_len/40) > 0 && $read_len/40 <4);
	
	my $soap_fwd_map_file = $fqName.".fwd.soap.out";
	my $cmd = "$soap2_exe -r 0 -v $maxMismatches -p $cpu -D $template_fwd -a $encodedFqName -o $soap_fwd_map_file 2>&1";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
	print "$cmd\n";
	my $soap_rev_map_file = $fqName.".rev.soap.out";
	my $cmd = "$soap2_exe -r 0 -v $maxMismatches -p $cpu -D $template_rev -a $encodedFqName -o $soap_rev_map_file 2>&1";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
	print "$cmd\n";
	unlink($encodedFqName);

	return ($soap_fwd_map_file, $soap_rev_map_file);
}

sub fastq2SOAPpe(){
	my $tmp1 = shift;
	my $tmp2 = shift;
	my @f = split /\//, $tmp1;
	my $fqName1 = pop(@f);
	@f = split /\//, $tmp2;
	my $fqName2 = pop(@f);
	my $encodedFqName1 = $fqName1 . ".encoded";
	my $encodedFqName2 = $fqName2 . ".encoded";
	encodeFastq($tmp1, $encodedFqName1);
	encodeFastq($tmp2, $encodedFqName2);

	my $maxMismatches=1;
	$maxMismatches = int($read_len/40) if(int($read_len/40) > 0 && $read_len/40 <4);
	
	my $soap_fwd_map_PE_file = $fqName1.".fwd.soap.PE.out";
	my $soap_fwd_map_SE_file = $fqName1.".fwd.soap.SE.out";
	my $cmd = "$soap2_exe -r 0 -v $maxMismatches -p $cpu -D $template_fwd -a $encodedFqName1 -b $encodedFqName2 -o $soap_fwd_map_PE_file -2 $soap_fwd_map_SE_file -m 1 -x 400 2>&1";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
	print "$cmd\n";
	my $soap_rev_map_PE_file = $fqName1.".rev.soap.PE.out";
	my $soap_rev_map_SE_file = $fqName1.".rev.soap.SE.out";
	my $cmd = "$soap2_exe -r 0 -v $maxMismatches -p $cpu -D $template_rev -a $encodedFqName1 -b $encodedFqName2 -o $soap_rev_map_PE_file -2 $soap_rev_map_SE_file -m 1 -x 400 2>&1";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
	print "$cmd\n";
	unlink($encodedFqName1);
	unlink($encodedFqName2);

	system("cat $soap_fwd_map_SE_file >> $soap_fwd_map_PE_file");
	system("cat $soap_rev_map_SE_file >> $soap_rev_map_PE_file");

	unlink($soap_fwd_map_SE_file);
	unlink($soap_rev_map_SE_file);

	return ($soap_fwd_map_PE_file, $soap_rev_map_PE_file);

}

sub makeCalls(){
        my $strand = shift;
        my $hash = shift;
        my $array = shift;
        my $cpg_h = shift;
        my $snp_h = shift;
        my %cpgTable = %$hash;
        my @fields = split /\t/, $array;
        my $pos = $fields[1];
        $pos = $chrSizes{$fields[0]}-$fields[1]+1 if($strand eq "C");
        my $index = $fields[0] . ":". $strand . "\t". $pos;
        if($cpgTable{$index}){  #this is any C site
                my %variantStat = pileupFields2variantStat(\@fields, 0);
                print $cpg_h $fields[0], "\t", $pos, "\t", $strand, "\t", $variantStat{"depth"}, "\t", $cpgTable{$index};
                foreach my $base (keys(%{$variantStat{"counts"}})){
                        print $cpg_h "\t$base\t", $variantStat{"counts"}->{$base};
                }
                print $cpg_h "\n";
        }
        if($fields[4] ne '.' && $fields[5] > 20 && $snp_file){  # if the SNP quality is higher than 20
                my %variantStat = pileupFields2variantStat(\@fields, 0);
                print $snp_h $fields[0], "\t", $pos, "\t", $variantStat{"refBase"}, "\t", $strand, "\t", $variantStat{"call"}, "\t", $variantStat{"snpQual"}, "\t", $variantStat{"depth"};
                foreach my $base (keys(%{$variantStat{"counts"}})){
                        print $snp_h "\t$base\t", $variantStat{"counts"}->{$base};
                }
                print $snp_h "\n";
        }
        return if($strand eq "W"); # skip the following if strand is forward.
        $pos--; #cpg pos Crick to Watson
        $index = $fields[0] . ":W\t". $pos;
        return if(!$cpgTable{$index});
        if($cpgTable{$index} eq "CG"){ #this is a CG site indexed on Watson
                my %variantStat = pileupFields2variantStat(\@fields, 0);
                print $cpg_h $fields[0], "\t", $pos, "\tC\t", $variantStat{"depth"}, "\t", $cpgTable{$index};
                foreach my $base (keys(%{$variantStat{"counts"}})){
                        print $cpg_h "\t$base\t", $variantStat{"counts"}->{$base};
                }
                print $cpg_h "\n";
        }
        return;
}


	
sub pileupFields2variantStat(){	
	my ($h_fields, $mask_methyl) = @_;
	my @fields = @{$h_fields};
	my $refBase = $fields[2];
	my $readBase = $fields[8];
	my $readQual = $fields[9];
	$readBase =~ s/\$//g;
	$readBase =~ s/\^//g;
	$readBase =~ s/F//g;
	$readBase =~ s/-[0-9]+[ACGTNacgtn]+/D/g;
	$readBase =~ s/-[0-9]+[ACGTNacgtn]+/I/g;
	$readBase =~ s/[\.\,]/$refBase/g;
	my $totalCounts=0;
	my %variantStat;
	while(my $base = chop($readBase)){
		my $baseQual = chop($readQual); 
		next if($base !~/[ATGC]/i || ord($baseQual)-$qual_base < 5); #ignore low-quality bases (phred score < 5)   
		$base=uc($base);
		$base =~ s/C/T/g if($mask_methyl); # mask methylation
		$variantStat{'counts'}->{$base}++;	
		$totalCounts++;		
	}
	$variantStat{'refBase'}=$refBase;
	$variantStat{'depth'}= $totalCounts;
	$variantStat{'snpQual'}= $fields[5];
	$variantStat{'call'}= $fields[3];
	return %variantStat;
}


sub encodeFastq(){
	my $fileName = shift;
	my $outFileName = shift;
	open(FQ, "$fileName") || die ("Error in opening file $fileName!");
	open(FQ_OUT, ">$outFileName") || die ("Error in opening file $outFileName!");
	while(my $line1 = <FQ>){
		chop($line1);
		my $line2 = <FQ>; 
		chop($line2);
		my $line3 = <FQ>;
		chop($line3);
		my $line4 = <FQ>;
		chop($line4);
		if($threep || $fivep){
			my $start = $fivep;
			my $total = length($line2) - $threep - $fivep;
			my $tmp = substr($line2, $start, $total);
			$line2=$tmp;
			$tmp = substr($line4, $start, $total);
			$line4=$tmp;
		}
		if(guess_strand($line2) eq "R"){			
			my $seq = revComp($line2);
			my $qual ="";
			while(my $qualBase = chop($line4)){
				$qual = $qual . $qualBase;
			}
			$line2=$seq;
			$line4=$qual;			
		}
		$line1 = $line1 . "|" . $line2;
		$line3 = $line3 . "|" . $line2;
		$line2 =~ s/C/T/g;
		print FQ_OUT "$line1\n$line2\n$line3\n$line4\n";
	}
	close(FQ);
	close(FQ_OUT);
}

sub revComp(){
	my $seq = shift;
	my $rcSeq='';
	for(my $i=0; $i<length($seq); $i++){
		$rcSeq = $rcTable{substr($seq,$i,1)} . $rcSeq;
	}
	return $rcSeq;
}

sub guess_strand(){
	my $seq = shift;
	my %baseCounts;
	$baseCounts{'A'}=0.001;
	$baseCounts{'T'}=0.001;
	$baseCounts{'G'}=0.001;
	$baseCounts{'C'}=0.001;
	while(my $base = chop($seq)){
		$baseCounts{$base}++;
	}
	if($baseCounts{'T'}/$baseCounts{'C'} >$baseCounts{'A'}/$baseCounts{'G'}) {
		return "F";
	}else{
		return "R";
	}
}
main();
