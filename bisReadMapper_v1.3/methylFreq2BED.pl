#!/usr/bin/perl -w
use strict;
# A perl script for converting a methylFreq file to a UCSC BED file. 
# USAGE: ./methylFreq2BED.pl sampleID [minDepth] [snp_file ] < sample_methylFreq.txt
#   sample_methylFreq.txt: the methylFreq file generated by bisReadMapper.pl
# Format of the BED file:
# Column 1: chromosome
# Column 2: chromosome_start
# Column 3: chromosome_end
# Column 4: methylation level [0-1]
# column 5: sequencing depth
# column 6: strand, please ignored since the reads from both strands were combined
# Column 7: chromosome_start
# Column 8: chromosome_end
# Column 9: pesudo color, red is methylated, and green is unmethylated
#
# Written by Kun Zhang (kzhang@bioeng.ucsd.edu), last modified 10/07/2010 
#

my $sampleID = $ARGV[0];
$sampleID = "Sample" if (!$sampleID);
my $minDepth = $ARGV[1];
$minDepth = 10 if(!$minDepth);
my $snp_file = $ARGV[2];

my $cg_file = $sampleID . ".cg-pos.BED.txt";
my $chg_file = $sampleID . ".chg-pos.BED.txt";
my $chh_file = $sampleID . ".chh-pos.BED.txt";
my $c_snp_file = $sampleID . ".cSNP-pos.BED.txt" if($snp_file);
my $chg_called = 0;
my $chh_called = 0;

my @palette=("0,240,0", "30,210,0","60,180,0","90,150,0","120,120,0","150,90,0","180,60,0","210,0,0");

my %methylTable;

sub main(){
	while(my $line = <STDIN>){
		chop($line);
		my @fields = split(/\t/, $line);		
		my $strand = $fields[2] eq 'W' ? '+' : '-';
		my %alleleCounts;
		my $CT_counts;
		$fields[1]-- if($strand eq '-');
		for(my $i=5; $i<scalar(@fields); $i+=2){
			$alleleCounts{$fields[$i]}=$fields[$i+1];			
			$CT_counts += $fields[$i+1] if($fields[$i]=~ /[CT]/);
		}
		next if(!$CT_counts || $CT_counts/$fields[3] < 0.9);
		my $index=$fields[0] . ":" . $fields[1];
		$alleleCounts{'C'} =0 if(!$alleleCounts{'C'});
		$methylTable{$index}->{'C'} +=  $alleleCounts{'C'} ;
		$methylTable{$index}->{'CT'} += $CT_counts;
		$methylTable{$index}->{"Type"} = $fields[4];
		$chg_called = 1 if($fields[4] eq 'CHG');
		$chh_called = 1 if($fields[4] eq 'CHH');
		$methylTable{$index}->{"Strand"} = $strand;
	}
	load_dbsnp($snp_file) if($snp_file);
	report_methylFreqBED();
}

sub report_methylFreqBED(){
	open(CG_OUT, ">$cg_file") || die("Error writing CG BED file\n");
	open(CHH_OUT, ">$chh_file") || die("Error writing CHH BED file\n") if($chh_called);
	open(CHG_OUT, ">$chg_file") || die("Error writing CHG BED file\n") if($chg_called);
	open(C_SNP_OUT, ">$c_snp_file") || die("Error writng C-SNP BED file\n") if($snp_file);
	print CG_OUT "track name=\"", $sampleID, "\" description=\"Methylation level at CGs\" visibility=2 useScore=1 itemRgb=\"On\"\n";
	print CHH_OUT "track name=\"", $sampleID, "\" description=\"Methylation level at CHH\" visibility=2 useScore=1 itemRgb=\"On\"\n" if($chh_called);
	print CHG_OUT "track name=\"", $sampleID, "\" description=\"Methylation level at CHG\" visibility=2 useScore=1 itemRgb=\"On\"\n" if($chg_called);
	print C_SNP_OUT "track name=\"", $sampleID, "\" description=\"Methylation level at C-SNP\" visibility=2 useScore=1 itemRgb=\"On\"\n" if($snp_file);
	foreach my $index(sort keys(%methylTable)){
		next if($methylTable{$index}->{'CT'}<$minDepth);
		my ($chr,$chr_pos) = split(/:/, $index);	
		my $methylLevel = sprintf("%4.3f", $methylTable{$index}->{'C'} ? $methylTable{$index}->{'C'}/$methylTable{$index}->{'CT'} : 0);
		if($methylTable{$index}->{"rs"}){
			print C_SNP_OUT "$chr\t", 
				$chr_pos-1, "\t", 
				$chr_pos, "\t",
				"$methylLevel\t",
				int($methylTable{$index}->{'CT'}), "\t", 
				$methylTable{$index}->{"Strand"}, "\t",
				$chr_pos-1, "\t",
				$chr_pos, "\t",
				$palette[int($methylLevel*8-0.0001)],"\n";
			next;
		}else{
			print CG_OUT "$chr\t", 
				$chr_pos-1, "\t", 
				$chr_pos, "\t",
				"$methylLevel\t",
				int($methylTable{$index}->{'CT'}), "\t+\t",
				$chr_pos-1, "\t",
				$chr_pos, "\t",
				$palette[int($methylLevel*8-0.0001)],"\n" if($methylTable{$index}->{"Type"} eq 'CG');
			print CHG_OUT "$chr\t", 
				$chr_pos-1, "\t", 
				$chr_pos, "\t",
				"$methylLevel\t",
				int($methylTable{$index}->{'CT'}), "\t", 
				$methylTable{$index}->{"Strand"}, "\t",
				$chr_pos-1, "\t",
				$chr_pos, "\t",
				$palette[int($methylLevel*8-0.0001)],"\n" if($methylTable{$index}->{"Type"} eq 'CHG');
			print CHH_OUT "$chr\t", 
				$chr_pos-1, "\t", 
				$chr_pos, "\t",
				"$methylLevel\t",
				int($methylTable{$index}->{'CT'}), "\t", 
				$methylTable{$index}->{"Strand"}, "\t",
				$chr_pos-1, "\t",
				$chr_pos, "\t",
				$palette[int($methylLevel*8-0.0001)],"\n" if($methylTable{$index}->{"Type"} eq 'CHH');

		}
		
	}
	close(CG_OUT);
	close(CHH_OUT);
	close(CHG_OUT);
	close(C_SNP_OUT);
}

sub load_dbsnp(){
	my $refSnpFile = shift;
	open(INFILE, "$refSnpFile") || die("Error in opening $refSnpFile!\n");
	while(my $line = <INFILE>){
		my @fields = split(/\t/, $line);
		my $index = $fields[1] . ":" . $fields[3];
		$fields[3]++;
		my $index_g = $fields[1] . ":" . $fields[3];
		if($methylTable{$index}){
			if($fields[6] eq $methylTable{$index}->{"Strand"}) {
				$methylTable{$index}->{"rs"} = $fields[4];
			}elsif($methylTable{$index}->{"Type"} eq 'CG'){
				$methylTable{$index}->{"rs"} = $fields[4];
			}		
		}elsif($methylTable{$index_g}){
			$methylTable{$index_g}->{"rs"} = $fields[4] if($methylTable{$index_g}->{"Type"} eq 'CG');
		}
	}
	close(INFILE);
}

main();
