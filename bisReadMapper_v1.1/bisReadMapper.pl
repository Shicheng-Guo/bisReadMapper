#!/usr/bin/perl
use strict;
use Switch;
# bisReadMapper.pl: a perl script to map single-end bisulfite sequencing reads and report the methylation levels and SNPs. 
# USAGE: ./bisReadMapper.pl params.txt > SAMPLE_NAME.log
# This version performs soft trimming and the old samtools pileup command.

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
my $qualtrim=0;

my $script_dir = `readlink -f $0`;
chomp($script_dir);
$script_dir =~ s/\/bisReadMapper.pl//g;

my $keep_bam = 1;
my $minDepth = 10;
my $rmdup = 0;
my $num_lines = $ARGV[1];
$num_lines = 10**8 if(!$num_lines);

my ($template_fwd, $template_fwd_fa, $template_rev, $template_rev_fa, $template_idx);

my ($soap2_exe, $samtools, $soap2sam);

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
my %chrFiles;
my %samList;
my %mappedReads;

sub main(){
	open(PARAMS, "$ARGV[0]") || die("No params file given or params file is unreadable\n");
	while(my $line = <PARAMS>){
		chomp($line);
		next if($line =~ /^#/);
		my ($arg, $val) = split("=", $line);
		switch($arg) {
			case "reads" { @reads = split(",", $val); print "Reads: ", join(",", @reads), "\n"; }
			case "length" { $read_len = int($val); print "Read length to use: $read_len\n"; }
			case "refDir" { $ref_dir = $val; print "Reference dir: $ref_dir\n"; }
			case "soapDir" { $soap_dir = $val; print "SOAP dir: $soap_dir\n"; }
			case "alignMode" { $align_mode = $val; print "Reads alignment mode: $align_mode\n"; }
			case "qualBase" { $qual_base = int($val); print "Base quality: $qual_base\n"; }
			case "numCPU" { $cpu = int($val); print "Number of processors for mapping: $cpu\n"; }
			case "samtoolsDir" { $samtools_dir = $val; print "Samtools dir: $samtools_dir\n"; }
			case "soap2sam" { $soap2sam = $val; print "Soap2sam path: $soap2sam\n"; }
			case "snp" { $snp_file = $val; print "dbSNP file: $snp_file\n"; }
			case "name" { $name = $val; print "Sample name: $name\n"; }
			case "allC" { $allC = isyes($val); print "Call all C? : $allC\n"; }
			case "trim3" { $threep = int($val); print "Bases to trim from 3': $threep\n"; }
			case "trim5" { $fivep = int($val); print "Bases to trim from 5': $fivep\n"; }
			case "bam" { $bam = isyes($val); print "Input is bam? : $bam\n"; }
			case "qualTrim" { $qualtrim = int($val); print "Perform quality trimming using cutoff : $qualtrim\n"; }
			#if bam, two reads file must be provided, reads from Watson first, then Crick strand.
			#assumes bam file have already been filtered and sorted.
			case "minDepth" { $minDepth = int($val); print "Minimum depth for BED: $minDepth\n"; }
			case "rmdup" { $rmdup = isyes($val); print "Remove duplicates (rmdup) ?: $rmdup\n"; }
		}
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
		die("Can't align by paired-end, only one reads file detected.\n");
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
	$template_idx = $template_fwd_fa. ".fai";
	
	open(GENOME_INDEX, "$template_idx") || die("Error opening chromosome sizes file for reading!\n");
	while(my $line = <GENOME_INDEX>){
		chomp($line);
		my @f = split "\t", $line;
		my $cur_chr = $f[0];
		$cur_chr =~ s/_Watson//;
		$cur_chr =~ s/_Crick//;
		$f[0] = $cur_chr;
		$samList{$cur_chr."_Watson.sam"} = 'W';
		$samList{$cur_chr."_Crick.sam"} = 'C';
		my $file = $cur_chr . "_Watson.sam";
		open(TMP_OUT, ">$file") || die("Error writing to $file\n");
		close(TMP_OUT);
		$file = $cur_chr . "_Crick.sam";
		open(TMP_OUT, ">$file") || die("Error writing to $file\n");
		close(TMP_OUT);
		$chrSizes{$cur_chr} = $f[1];
		# find the chr position files:
		for my $val (@refList){
			next if($val != m/positions.txt/);
			if($val =~ m/$cur_chr.cpositions.txt/){
				push(@{$chrFiles{$cur_chr}},$val);
				print $cur_chr, "\t", $val, "\n";
				open(CPG_OUT, ">$name.$val.methylFreq") || die("Error writing to $name.$val.methylFreq.\n");
				close(CPG_OUT);
			}
		}
	}
	close(GENOME_INDEX);

        my $snpcall_file = $name . ".snp";
        open( my $snp_h, ">$snpcall_file") || die("Error writing to snp file, $snpcall_file \n");

	if($bam == 0){
		my ($soap_fwd_map_file, $soap_rev_map_file) = (0,0);
		if($align_mode eq "P"){
			for(my $i = 0; $i < scalar(@reads); $i+2){
				die("Uneven pairs of reads files given.\n") if(!$reads[$i+1]);
				($soap_fwd_map_file, $soap_rev_map_file) = fastq2SOAPpe($reads[$i], $reads[$i+1]);
				soap2sam($soap_fwd_map_file, $soap_rev_map_file);
			}
		}else{
			for(my $i = 0; $i < scalar(@reads); $i++){
				($soap_fwd_map_file, $soap_rev_map_file) = fastq2SOAPse($reads[$i]);
				soap2sam($soap_fwd_map_file, $soap_rev_map_file);
			}
		}
		print "Finished mapping reads.\n";
		my $sum_mapped = 0;
		my $sum_rmdupe = 0;
		foreach my $sam_file (keys %samList){
			my $processed_bam = sort_rmdup($sam_file);
			unlink($sam_file);
			extract($processed_bam, $samList{$sam_file}, $snp_h);
			$sum_mapped+= $mappedReads{$sam_file}->{"mapped"};
			$sum_rmdupe+= $mappedReads{$sam_file}->{"rmdup"};
			next if($mappedReads{$sam_file}->{"mapped"} == 0);
			print "Proportion of clonal reads removed: ", 
				sprintf("%4.3f", 1 - ($mappedReads{$sam_file}->{"rmdup"}/$mappedReads{$sam_file}->{"mapped"})), "\n";
		}
		print "++++++++++++++++++++++++++++++++++\n";
		print "Total proportion of clonal reads removed: ", sprintf("%4.3f", 1 - ($sum_rmdupe/$sum_mapped)), "\n";
		print "++++++++++++++++++++++++++++++++++\n";
	}else{
		extract($reads[0], 'W', $snp_h);
		extract($reads[1], 'C', $snp_h);
	}

	close($snp_h);
	undef %chrSizes;
	undef %chrFiles;
	undef %samList;
	undef %mappedReads;

	#filter SNPs given dbSNP file:
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
		
	#print "Calculating bisulfite data quality by correlating the forward and reverse methylation frequencies for CpGs at $minDepth minimum depth \n";
	#$cmd = "less $name.*.methylFreq | $script_dir/frMethylCorr.pl $minDepth";
	#print $cmd, "\n";
	#system($cmd) == 0 or die "system problem (exit $?): $!\n";
	
}

sub sort_rmdup(){
        my $sam_file = shift;
        my $sorted_bam = $sam_file . ".sorted";
        $sorted_bam =~ s/sam/bam/;
        my $template = $template_fwd_fa;
        $template = $template_rev_fa if($sam_file =~ m/Crick/);

	#begin modified codes
	#if(!$rmdup){
	#        my $cmd = "$samtools view -ubST $template $sam_file | $samtools sort - $sorted_bam";
        #	print $cmd, "\n";
	#        system($cmd) == 0 or die "system problem (exit $?): $!\n";
	#        $sorted_bam = $sorted_bam . ".bam";
	#	return $sorted_bam;
	#}

	#begin my own rmdup code:
	#my $sorted_sam = $sam_file . ".sorted";
	#my $cmd = "sort -k4,4n $sam_file > $sorted_sam";
	#print $cmd, "\n";
	#system($cmd) == 0 or die "system problem (exit $?): $!\n";
	#my $rmdup_sam = "rmdup." . $sam_file;
	#open(SORTED, "$sorted_sam") || die("Error reading $sorted_sam\n");
	#open(RMDUP, ">$rmdup_sam") || die("Error writing to $rmdup_sam\n");
	#my $last_line = <SORTED>;
	#my @f = split "\t", $last_line;
	#my $last_pos = $f[3];
	#while(my $line = <SORTED>){
	#	@f = split "\t", $line;
	#	if($f[3] == $last_pos){
	#		$last_line = $line;
	#		$last_pos = $f[3];
	#	}else{
	#		#previous line is different from current line, print
	#		print RMDUP $last_line;
	#		$mappedReads{ $sam_file }->{"rmdup"}++;
	#		$last_line = $line;
	#		$last_pos = $f[3];
	#	}
	#}
	#print RMDUP $last_line;
	#$mappedReads{ $sam_file }->{"rmdup"}++;
	#close(SORTED);
	#unlink($sorted_sam);
	#unlink($sam_file);
        #$cmd = "$samtools view -ubST $template $rmdup_sam | $samtools sort - $sorted_bam";
        #print $cmd, "\n";
        #system($cmd) == 0 or die "system problem (exit $?): $!\n";
        #$sorted_bam = $sorted_bam . ".bam";
	#unlink($rmdup_sam);
	#return $sorted_bam;
	
        my $cmd = "$samtools view -ubST $template $sam_file | $samtools sort - $sorted_bam";
        print $cmd, "\n";
        system($cmd) == 0 or die "system problem (exit $?): $!\n";
	$sorted_bam = $sorted_bam . ".bam";
        my $sorted_rmdup_bam = "rmdup." . $sorted_bam;
        $cmd = "$samtools rmdup -S $sorted_bam $sorted_rmdup_bam";
        print $cmd, "\n";
        system($cmd) == 0 or die "system problem (exit $?): $!\n";
	unlink($sorted_bam);
	$cmd = "$samtools flagstat $sorted_rmdup_bam";
	my @tmp = split /\n/, `$cmd`;
	my @tmp2 = split " ", $tmp[0];
	$mappedReads{ $sam_file }->{"rmdup"} = int($tmp2[0]);
        return ($sorted_rmdup_bam);
}

sub soap2sam(){
	my ($fwd_map, $rev_map) = @_;
	#merge the two soap mapped files
	open(SOAP_REV, "$rev_map") || die("Error in opening $rev_map.");
	open(SOAP_FWD, ">>$fwd_map") || die("Error in opening $fwd_map.");
	while(my $line = <SOAP_REV>){
		chomp($line);
		print SOAP_FWD $line, "\n";
	}
	close(SOAP_REV);
	close(SOAP_FWD);	
	unlink($rev_map);
	
	my $fqName = $reads[0];
	#identify unique reads and save them in two files based on the templates.
	my $soap_combined_sorted_map_file = $fqName.".combined.sorted.soap.out";
	my $cmd = "sort -k 1,1 < $fwd_map > $soap_combined_sorted_map_file";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
	unlink($fwd_map);
	
	open(SOAP_OUT, "$soap_combined_sorted_map_file") || die("Error in opening $soap_combined_sorted_map_file.");

	my %file_handles;
	foreach my $val (keys %samList){
		open( $file_handles{$val} , ">$val.tmp" ) || die ("Error writing to file $val.tmp\n");
	}	

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
			print { $file_handles{ $last_fields[7] . ".sam"} } join("\t", @last_fields), "\n";
			$mappedReads{$last_fields[7] . ".sam"}->{"mapped"}++;
			$last_line = $line;
			@last_fields = @fields;
		}
	}
	#print the last line
	if($last_line){
		my ($id, $orig_seq) = split (/\|/, $last_fields[0]);
		$last_fields[0] = $id;
		$last_fields[1] = $orig_seq;
		print { $file_handles{$last_fields[7] . ".sam"} } join("\t", @last_fields), "\n";
		$mappedReads{$last_fields[7]. ".sam"}->{"mapped"}++;
	}
	close(SOAP_OUT);
	
	foreach my $val (keys %samList){
		close( $file_handles{$val} );
		my $cmd = "$soap2sam $val.tmp >> $val";
		print $cmd, "\n";
		system($cmd) == 0 or die "system problem (exit $?): $!\n";
		unlink($val . ".tmp");
	}
	unlink($soap_combined_sorted_map_file);
}

sub extract(){
        my ($sorted_bam, $str, $snp_h) = @_;

        my $cmd = "$samtools index $sorted_bam";
        system($cmd) == 0 or die "system problem (exit $?): $!\n";
        my $pileup = $sorted_bam . ".pileup";
	
	if($str eq 'W'){
	        $cmd = "$samtools pileup -cf $template_fwd_fa $sorted_bam > $pileup";
        	print $cmd, "\n";
	        system($cmd) == 0 or die "system problem (exit $?): $!\n";
	}else{
	        $cmd = "$samtools pileup -cf $template_rev_fa $sorted_bam > $pileup";
        	print $cmd, "\n";
	        system($cmd) == 0 or die "system problem (exit $?): $!\n";
	}


	open(FWD_FILE, "$pileup") || die("Error opening $pileup\n");
	my @cur_pos;
	my $cnt = 0;
	while(my $line = <FWD_FILE>){
		next if($line =~ m/##/);
		chomp($line);
		$cur_pos[$cnt] = $line;
		if($cnt == $num_lines){
			print "Processing $sorted_bam data\n";
			processCur(\@cur_pos, $cnt, $snp_h, $str);
			$cnt = 0;
		}else{
			$cnt++;
		}
	}
	close(FWD_FILE);
	#process the last positions
	if($cnt != 0){
		print "Processing last $sorted_bam data\n";
		processCur(\@cur_pos, $cnt-1, $snp_h, $str);
	}	

        unlink($sorted_bam . ".bai");
        unlink($pileup);

        return;
}

sub processCur(){
	my ($arrayAddy, $cnt, $snp_h, $str) = @_;
	my @array = @{$arrayAddy};
	my %posTable;
	my %chr_list;
	my $max_pos = 0;
	print "- $cnt lines\n";
	for(my $i = 0; $i<=$cnt; $i++){
		my $line = $array[$i];
		my @fields = split(/\t/, $line);
		my $chr = $fields[0];
		$chr =~ s/_Watson//;
		$chr =~ s/_Crick//;
		my $fwd_pos = $fields[1];
		$fwd_pos = $chrSizes{$chr}-$fields[1]+1 if($str eq 'C');
			
		#saves chromosomes to chr_list
		die("$chr\n$line\n") if(!$chrFiles{$chr});
		foreach my $file (@{$chrFiles{$chr}}){
			$chr_list{$file}=1;
		}
		#saves to posTable
		$posTable{$chr.":".$str.":".$fwd_pos} = $line;
		#set max position
		$max_pos = $fwd_pos if($fwd_pos > $max_pos);
	}

	print "~ Searching in ", join(",", keys %chr_list), "\n";
	foreach my $chr_file (keys %chr_list){
		my $methylFreq_file = $name . "." . $chr_file . ".methylFreq";
		open( my $cpg_h, ">>$methylFreq_file" ) || die("Error writing to $methylFreq_file");
		open(C_POS, "$ref_dir/$chr_file") || die("Error opening file $chr_file\n");
		while(my $line = <C_POS>){
			chomp($line);
			my @fields = split(/\t/, $line);
			my $pos = $fields[1];
			my $context = $fields[2];
			next if($context ne 'CG' && $allC == 0);
			my $chr_str = $fields[0];
			#don't continue if reached maximum position to search
			last if($pos > $max_pos);
			#now look up the position
			#if the position is a CG and the current positions are from Crick strand, only look up ($pos+1).
			if($context eq 'CG' && $str eq 'C'){
				$chr_str =~ s/:W/:C/;
				$pos = $pos++;
				if($posTable{$chr_str.":".$pos}){
					$fields[0] =~ s/:W//;
					$fields[0] =~ s/:C//;
					my @pileup_dat = split(/\t/, $posTable{$chr_str.":".$pos});
					next if($pileup_dat[0] eq "NA");
					$posTable{$chr_str.":".$pos} = "NA";
					my %variantStat = pileupFields2variantStat(\@pileup_dat, 0);
					next if($variantStat{"depth"} == 0);
					print $cpg_h $fields[0], "\t$pos\t$str\t", $variantStat{"depth"}, "\t", $context;
					foreach my $base (keys(%{$variantStat{"counts"}})){
						print $cpg_h "\t$base\t", $variantStat{"counts"}->{$base};
					}
					print $cpg_h "\n";
				}
			}else{
				#this is a CG/CHG/CHH position
				if($posTable{$chr_str.":".$pos}){
					$fields[0] =~ s/:W//;
					$fields[0] =~ s/:C//;
					my @pileup_dat = split(/\t/, $posTable{$chr_str.":".$pos});
					next if($pileup_dat[0] eq "NA");
					$posTable{$chr_str. ":". $pos} = "NA";
					my %variantStat = pileupFields2variantStat(\@pileup_dat, 0);
					next if($variantStat{"depth"} == 0);
					print $cpg_h $fields[0], "\t$pos\t$str\t", $variantStat{"depth"}, "\t", $context;
					foreach my $base (keys(%{$variantStat{"counts"}})){
						print $cpg_h "\t$base\t", $variantStat{"counts"}->{$base};
					}
					print $cpg_h "\n";
				}
			}
		}
		close(C_POS);
		close($cpg_h);
	}

	foreach my $val(keys %posTable){
		next if($posTable{$val} eq "NA");
		my @fields = split "\t", $posTable{$val};
		my $chr = $fields[0];
		$chr =~ s/_Watson//;
		$chr =~ s/_Crick//;
                my $fwd_pos = $fields[1];
                $fwd_pos = $chrSizes{$chr}-$fields[1]+1 if($str eq 'C');
		if($fields[4] ne '.' && $fields[5] > 20){  # if the SNP quality is higher than 20
                        my %variantStat = pileupFields2variantStat(\@fields, 0);
                        print $snp_h $chr, "\t$fwd_pos\t", $variantStat{"refBase"}, "\t$str\t",
                                $variantStat{"call"}, "\t", $variantStat{"snpQual"}, "\t", $variantStat{"depth"};
                        foreach my $base (keys(%{$variantStat{"counts"}})){
                                print $snp_h "\t$base\t", $variantStat{"counts"}->{$base};
                        }
                        print $snp_h "\n";
                }

	}
	undef %posTable;
	undef @array;
	undef %chr_list;
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
	$maxMismatches = int($read_len/30) if(int($read_len/30) > 0 && $read_len/30 <4);
	
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
		chomp($line1);
		my $line2 = <FQ>; 
		chomp($line2);
		my $line3 = <FQ>;
		chomp($line3);
		my $line4 = <FQ>;
		chomp($line4);
		if($threep || $fivep){
			my $start = $fivep;
			my $total = length($line2) - $threep - $fivep;
			$line2 = substr($line2, $start, $total);
			#$line2 = substr($line2, 0, 36);
			$line4 = substr($line4, $start, $total);
			#$line4 = substr($line4, 0, 36);
		}elsif($qualtrim){
		#http://wiki.bioinformatics.ucdavis.edu/index.php/TrimBWAstyle.pl
			my @quals = split "", $line4;
			my $pos = scalar(@quals);
			my $maxPos = $pos;
			my $sum_bad_qual = 0;
			my $maxSum = 0;
			while( $pos > 0 && $sum_bad_qual >= 0){
				$sum_bad_qual += $qualtrim - (ord($quals[$pos-1])-$qual_base);
				if($sum_bad_qual > $maxSum){
					$maxSum = $sum_bad_qual;
					$maxPos = $pos;
				}
				$pos--;
			}
			if($pos==0 || $maxPos < 30) { $line2 = "N\n"; $line4 = "B\n"; }
			else{
				$line2 = substr($line2,0,$maxPos);
				$line4 = substr($line4,0,$maxPos);
			}
		}
		#$line1 = $line1 . "|" . $line2;
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
	if($baseCounts{'T'}/$baseCounts{'C'} > $baseCounts{'A'}/$baseCounts{'G'}) {
		return "F";
	}else{
		return "R";
	}
}

sub isyes(){
	my $val = shift;
	if($val =~ m/^y/i) { return 1; } else{ return 0; }
}
main();
