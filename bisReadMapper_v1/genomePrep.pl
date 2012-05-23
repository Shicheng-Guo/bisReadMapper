#!/usr/bin/perl -w
# Usage 1: genomePrep.pl genome.fa context=[cg/all] convert=[yes/no]
# convert: yes = convert genome, no = don't convert genome
# context: CG = CG only, ALL = all C context

use strict;

my $genome_path = $ARGV[0];
die("Missing genome name!") if(!$ARGV[0]);
my @tmp = split /\//, $genome_path;
my $genome = pop(@tmp);
undef @tmp;

my $variant_file = $ARGV[1];

my $all = 0;
my $convert = 1;

for(my $i = 1; $i<scalar(@ARGV); $i++){
	my $val = $ARGV[$i];
	if($val =~ m/convert/){ $val=~s/convert=//g; $convert = 0 if($val =~ m/no/i); }
	if($val =~ m/context/){ $val=~s/context=//g; $all = 1 if($val =~ m/all/i); }
}

my %rcTable;
$rcTable{'A'}='T';
$rcTable{'T'}='A';
$rcTable{'G'}='C';
$rcTable{'C'}='G';

my %variant;

if(open(IN, "$variant_file")){
	$genome = "var." . $genome;
	while(my $line = <IN>){
		chop($line);
		my @f = split /\t/, $line;
		next if($f[10] ne "ref" && $f[16] ne "ref");
		$variant{$f[2]}->{$f[4]} = $f[11] if($f[10] eq "snp");
		$variant{$f[2]}->{$f[4]} = $f[17] if($f[16] eq "snp");
		#print $f[3], "\t", $f[5], "\t",  $f[10], "\t", $f[16], "\t", $f[11], "\t", $f[17], "\n";
	}
	close(IN);
}

my $bisFwd = $genome . ".bis.fwd";
my $bisRev = $genome . ".bis.rev";

if($convert){
        open(FWD_OUT, ">$bisFwd") || die("Error writing bis fwd file");
        close(FWD_OUT);
        open(REV_OUT, ">$bisRev") || die("Error writing bis rev file");
        close(REV_OUT);
}

if(open(IN, "$genome_path")){
	my $line = <IN>;
	chop($line);
	while($line =~ m/>/){
		$line =~ s/>//g;	
		my $name = $line;
		#die if($name ne "chr10"); ###remember to remove!!
		my $rev_seq = "\n";
		my $fwd_seq = "NA";
		my @rcArray;
		my $pos = 0;
		my ($first, $second, $third) = ('N', 'N', 'N');
		open(CH_OUT, ">$genome.$name.cpositions.txt") || die("Error writing c positions file");

		if($convert){
			open(FWD_OUT, ">>$bisFwd") || die("Error writing bis fwd file");
			print FWD_OUT ">", $name, "\n";
		}

		while(my $seq = <IN>){	
			chop($seq);
			$seq =~ tr/ //;
			next if(!$seq);
			if($seq =~ m/>/){
				$line = $seq;
				last;
			}
			$seq = uc($seq);
			my @base = split("", $seq);
			for(my $i = 0; $i < scalar(@base); $i++){
				$pos++;
				my $val = $base[$i];
				if($variant{$name}->{$pos}){
					$val = $variant{$name}->{$pos};
					$base[$i] = $val;
				}
				my $rcval = $val;
				$rcval = $rcTable{$val} if($rcTable{$val});
				$rev_seq = $rcval . $rev_seq;
				if($fwd_seq eq "NA"){
					$fwd_seq = $val;
				}else{
					$fwd_seq = $fwd_seq . $val;
				}
				if($pos > 3){
					my $c_pos = $pos - 3;
					print CH_OUT $name, ":W\t", $c_pos, "\tCG\n" if($first eq 'C' and $second eq 'G');
					if($all){
						if($first eq 'C' and $second ne 'G'){
							print CH_OUT $name, ":W\t", $c_pos, "\tCHG\n", $name, ":C\t", $c_pos+2, "\tCHG\n" 
								if($third eq 'G');
							print CH_OUT $name, ":W\t", $c_pos, "\tCHH\n" 
								if($third ne 'G');
						}
						print CH_OUT $name, ":C\t", $c_pos+2, "\tCHH\n" 
							if($first ne 'C' and $second ne 'C' and $third eq 'G');
					}
				}
				($first, $second, $third) = ($second, $third, $val);
				if($convert && $pos%200 == 0){
					$fwd_seq =~ s/C/T/g;
					$rev_seq =~ s/C/T/g;
					print FWD_OUT $fwd_seq, "\n";
					unshift(@rcArray, $rev_seq);
					$rev_seq = "\n";
					$fwd_seq = "NA";
				}
			}
			undef @base;
		}
		close(CH_OUT);

		if($convert){
			if($pos%200 != 0){
				$fwd_seq =~ s/C/T/g;
				$rev_seq =~ s/C/T/g;
				print FWD_OUT $fwd_seq, "\n";
				close(FWD_OUT);
				push(@rcArray, $rev_seq);
			}
			open(REV_OUT, ">>$bisRev") || die("Error writing bis rev file");
			print REV_OUT ">", $name, "\n";
			print REV_OUT @rcArray;
			close(REV_OUT);
		}
		undef @rcArray;
	}
	close(IN);
}
