#!/usr/bin/perl -w
# bisSnpFilter.pl: a perl script to combine raw SNP calls from both strands of DNA 
#                   and report confident SNP calls from bisulfite sequencing data. 
# USAGE: bisSnpFilter.pl  raw_snp_file dbSnp_file > filtered_snp_file
#

use strict;
my $snp132_file = $ARGV[1];

my %rcTable;
$rcTable{"A"}="T";
$rcTable{"T"}="A";
$rcTable{"G"}="C";
$rcTable{"C"}="G";
$rcTable{"N"}="N";
$rcTable{"R"}="Y";
$rcTable{"Y"}="R";
$rcTable{"M"}="K";
$rcTable{"K"}="M";
$rcTable{"S"}="S";
$rcTable{"W"}="W";
my %three2one;
$three2one{"A/G"}="R";
$three2one{"C/T"}="Y";
$three2one{"A/C"}="M";
$three2one{"G/T"}="K";
$three2one{"C/G"}="S";
$three2one{"A/T"}="W";

		my %callRuleW;              #	alleles	fwd	genotype
		$callRuleW{"R"}->{"R"}="R";	#	A/G(R)	A/G(R)	A/G(R)
		$callRuleW{"R"}->{"A"}="A";	#	A/G(R)	A		A
		$callRuleW{"R"}->{"G"}="G";	#	A/G(R)	G		G
		$callRuleW{"W"}->{"W"}="W";	#	A/T(W)	A/T(W)	A/T(W)
		$callRuleW{"W"}->{"A"}="A";	#	A/T(W)	A		A
		$callRuleW{"W"}->{"T"}="T";	#	A/T(W)	T		T
		$callRuleW{"M"}->{"W"}="M";	#	A/C(M)	A/T(W)	A/C(M)
		$callRuleW{"M"}->{"A"}="A";	#	A/C(M)	A		A
		$callRuleW{"M"}->{"T"}="C";	#	A/C(M)	T		C
		$callRuleW{"M"}->{"C"}="C";	#	A/C(M)	C		C
		$callRuleW{"K"}->{"K"}="K";	#	G/T(K)	G/T(K)	G/T(K)
		$callRuleW{"K"}->{"G"}="G";	#	G/T(K)	G		G
		$callRuleW{"K"}->{"T"}="T";	#	G/T(K)	T		T
		$callRuleW{"S"}->{"K"}="S";	#	C/G(S)	G/T(K)	C/G(S)
		$callRuleW{"S"}->{"G"}="G";	#	C/G(S)	G		G
		$callRuleW{"S"}->{"C"}="C";	#	C/G(S)	C		C
		$callRuleW{"S"}->{"T"}="C";	#	C/G(S)	T		C
		$callRuleW{"Y"}->{"Y"}="Y/C/T?";	#	C/T(Y)	C/T(Y)	Y/C/T?
		$callRuleW{"Y"}->{"T"}="Y/C/T?";	#	C/T(Y)	T		Y/C/T?
		$callRuleW{"Y"}->{"C"}="C";	#	C/T(Y)	C		C
		$callRuleW{"Y"}->{"R"}="R";	#	C/T(Y)	A/G(R) 	A/G(R)          SNP alleles on the reverse strand??
		$callRuleW{"M"}->{"Y"}="C";	#	A/C(M)	C/T(Y)	C/C(C)
		$callRuleW{"M"}->{"A"}="A";	#	A/C(M)	A		A
		$callRuleW{"M"}->{"T"}="C";	#	A/C(M)	T		C
		$callRuleW{"M"}->{"C"}="C";	#	A/C(M)	C		C
		$callRuleW{"M"}->{"M"}="M";	#	A/C(M)	A/C(M)	A/C(M)
		$callRuleW{"S"}->{"T"}="C";	#	C/G(S)	T		C
		$callRuleW{"S"}->{"K"}="S";	#	C/G(S)	T/G(K)	C/G(S)
		$callRuleW{"S"}->{"G"}="G";	#	C/G(S)	G		G

		my %callRuleC;              #   Ref     rev     genotype
		$callRuleC{"M"}->{"G"}="C";	#	A/C(M)	G		C
		$callRuleC{"M"}->{"T"}="A";	#	A/C(M)	T		A
		$callRuleC{"M"}->{"K"}="M";	#	A/C(M)	T/G(K)	A/C(M)
		$callRuleC{"M"}->{"W"}="M";	#	A/C(M)	A/T(W)	A/C(M)
		$callRuleC{"R"}->{"C"}="G";	#	A/G(R)	C		G
		$callRuleC{"R"}->{"Y"}="R/A/G?";	#	A/G(R)	C/T(Y)	R/A/G?
		$callRuleC{"R"}->{"T"}="R/A/G?";	#	A/G(R)	T		R/A/G?
		$callRuleC{"W"}->{"A"}="T";	#	A/T(W)	A		T
		$callRuleC{"W"}->{"W"}="W";	#	A/T(W)	A/T(W)	A/T(W)
		$callRuleC{"W"}->{"T"}="A";	#	A/T(W)	T		A
		$callRuleC{"S"}->{"C"}="G";	#	C/G(S)	C		G
		$callRuleC{"S"}->{"G"}="C";	#	C/G(S)	G		C
		$callRuleC{"S"}->{"K"}="S";	#	C/G(S)	G/T(K)	C/G(S)
		$callRuleC{"S"}->{"T"}="G";	#	C/G(S)	T		G
		$callRuleC{"Y"}->{"A"}="T";	#	C/T(Y)	A		T
		$callRuleC{"Y"}->{"G"}="C";	#	C/T(Y)	G		C
		$callRuleC{"Y"}->{"R"}="Y";	#	C/T(Y)	G/A(R)	C/T(Y)
		$callRuleC{"K"}->{"A"}="T";	#	G/T(K)	A		T
		$callRuleC{"K"}->{"W"}="K";	#	G/T(K)	A/T(W)	G/T(K)
		$callRuleC{"K"}->{"T"}="G";	#	G/T(K)	C		G
		$callRuleC{"K"}->{"Y"}="G";	#	G/T(K)	C/T(Y)	G
		$callRuleC{"K"}->{"C"}="G";	#	G/T(K)	T		G
		
		my %callRuleDS;		            #	fwd	    rev	    genotype
		$callRuleDS{"A"}->{"W"}="W";	#	A		A/T(W)	A/T(W)
		$callRuleDS{"A"}->{"Y"}="R";	#	A		C/T(Y)	A/G(R)
		$callRuleDS{"A"}->{"T"}="A";	#	A		T		A/A(A)
		$callRuleDS{"A"}->{"K"}="M";	#	A		T/G(K)	A/C(M)
		$callRuleDS{"R"}->{"Y"}="R";	#	A/G(R)	C/T(Y)	A/G(R)
		$callRuleDS{"R"}->{"T"}="R";	#	A/G(R)	T		A/G(R)
		$callRuleDS{"W"}->{"A"}="W";	#	A/T(W)	A		A/T(W)
		$callRuleDS{"W"}->{"W"}="W";	#	A/T(W)	A/T(W)	A/T(W)
		$callRuleDS{"W"}->{"G"}="M";	#	A/T(W)	G		A/C(M)
		$callRuleDS{"W"}->{"K"}="M";	#	A/T(W)	G/T(K)	A/C(M)
		$callRuleDS{"W"}->{"T"}="W/M?";	#	A/T(W)	T		A/T(W)?
		$callRuleDS{"W"}->{"T"}="W/M?";	#	A/T(W)	T		A/C(M)?
		$callRuleDS{"C"}->{"G"}="C";	#	C		G		C/C(C)
		$callRuleDS{"C"}->{"R"}="Y";	#	C		G/A(R)	C/T(Y)
		$callRuleDS{"C"}->{"K"}="M";	#	C		T/G(K)	A/C(M)
		$callRuleDS{"Y"}->{"R"}="Y";	#	C/T(Y)	G/A(R)	C/T(Y)
		$callRuleDS{"Y"}->{"G"}="Y";	#	C/T(Y)	G    	C/T(Y)
		$callRuleDS{"Y"}->{"A"}="Y";	#	C/T(Y)	A    	C/T(Y)
		$callRuleDS{"G"}->{"W"}="K";	#	G		A/T(W)	G/T(K)
		$callRuleDS{"G"}->{"C"}="G";	#	G		C		G/G(G)
		$callRuleDS{"G"}->{"Y"}="R/G?";	#	G		C/T(Y)	A/G(R)?
		$callRuleDS{"G"}->{"Y"}="R/G?";	#	G		C/T(Y)	G?
		$callRuleDS{"G"}->{"K"}="S";	#	G		G/T(K)	C/G(S)
		$callRuleDS{"G"}->{"T"}="G";	#	G		T		G/G(G)
		$callRuleDS{"K"}->{"A"}="K";	#	G/T(K)	A		G/T(K)
		$callRuleDS{"K"}->{"W"}="K";	#	G/T(K)	A/T(W)	G/T(K)
		$callRuleDS{"K"}->{"G"}="S";	#	G/T(K)	G		C/G(S)
		$callRuleDS{"K"}->{"K"}="S";	#	G/T(K)	G/T(K)	C/G(S)
		$callRuleDS{"K"}->{"T"}="K/S?";	#	G/T(K)	T		G/T(K)?
		$callRuleDS{"K"}->{"T"}="K/S?";	#	G/T(K)	T		C/G(S)?
		$callRuleDS{"T"}->{"A"}="T";	#	T		A		T/T(T)
		$callRuleDS{"T"}->{"W"}="W/K?";	#	T		A/T(W)	A/T(W)?
		$callRuleDS{"T"}->{"W"}="W/K?";	#	T		A/T(W)	G/T(K)?
		$callRuleDS{"T"}->{"G"}="C";	#	T		G		C/C(C)
		$callRuleDS{"T"}->{"R"}="Y";	#	T		G/A(R)	C/T(Y)
		$callRuleDS{"T"}->{"K"}="S/M?";	#	T		G/T(K)	C/G(S)?
		$callRuleDS{"T"}->{"K"}="S/M?";	#	T		G/T(K)	A/C(M)?
		
		my %callRuleDsRef;                     	#	fwd		rev		Ref			genotype
		$callRuleDsRef{"W"}->{"W"}->{"W"}="W";	#	A/T(W)	A/T(W)	W			W
		$callRuleDsRef{"W"}->{"T"}->{"W"}="W";	#	A/T(W)	T		W			W
		$callRuleDsRef{"W"}->{"T"}->{"M"}="M";	#	A/T(W)	T		M			M
		$callRuleDsRef{"W"}->{"T"}->{"M"}="M";	#	A		T		M			A
		$callRuleDsRef{"T"}->{"K"}->{"S"}="S";	#	T		G/T(K)	S			S
		$callRuleDsRef{"T"}->{"K"}->{"M"}="M";	#	T		G/T(K)	M			M
		$callRuleDsRef{"K"}->{"T"}->{"K"}="K";	#	G/T(K)	T		K			K
		$callRuleDsRef{"K"}->{"M"}->{"K"}="K";	#	G/T(K)	A/C(M)	G/T(K)		G/T(K)   
		$callRuleDsRef{"K"}->{"T"}->{"S"}="S";	#	G/T(K)	T		S			S
		$callRuleDsRef{"K"}->{"G"}->{"S"}="S";	#	G/T(K)	G		S			S
		$callRuleDsRef{"K"}->{"K"}->{"S"}="S";	#	G/T(K)	G/T(K)	S			S
		$callRuleDsRef{"K"}->{"S"}->{"S"}="S";	#	G/T(K)	C/G(S)	S			S
		$callRuleDsRef{"K"}->{"C"}->{"S"}="S";	#	G/T(K)	C		S			S
		$callRuleDsRef{"K"}->{"W"}->{"K"}="K";	#	G/T(K)	A/T(W)	G/T(K)		G/T(K)
		$callRuleDsRef{"T"}->{"W"}->{"W"}="W";	#	T		A/T(W)	W			W
		$callRuleDsRef{"T"}->{"W"}->{"K"}="K";	#	T		A/T(W)	K			K
		$callRuleDsRef{"T"}->{"A"}->{"K"}="T";	#	T		A		G/T(K)		T
		$callRuleDsRef{"R"}->{"Y"}->{"R"}="R";	#	A/G(R)	C/T(Y)	A/G(R)		A/G(R)
		$callRuleDsRef{"A"}->{"Y"}->{"R"}="R";	#	A    	C/T(Y)	A/G(R)		A/G(R)
		$callRuleDsRef{"R"}->{"C"}->{"R"}="R";	#	A/G(R)	C		A/G(R)		A/G(R)
		$callRuleDsRef{"G"}->{"Y"}->{"R"}="R";	#	G		C/T(Y)	A/G(R)		A/G(R)
		$callRuleDsRef{"G"}->{"Y"}->{"S"}="G";	#	G		C/T(Y)	S			G
		$callRuleDsRef{"G"}->{"Y"}->{"K"}="G";	#	G		C/T(Y)	K			G
		$callRuleDsRef{"G"}->{"Y"}->{"R"}="R/G?";#	G		C/T(Y)	A/G(R)		G/R?
		$callRuleDsRef{"G"}->{"C"}->{"R"}="G";	#	G		C		A/G(R)		G
		$callRuleDsRef{"G"}->{"C"}->{"S"}="G";	#	G		C		C/G(S)		G
		$callRuleDsRef{"G"}->{"C"}->{"K"}="G";	#	G		C		G/T(K)		G
		$callRuleDsRef{"G"}->{"T"}->{"K"}="G";	#	G		T		G/T(K)		G
		$callRuleDsRef{"G"}->{"T"}->{"S"}="G";	#	G		T		C/G(S)		G
		$callRuleDsRef{"Y"}->{"R"}->{"Y"}="Y";	#	C/T(Y)	A/G(R)	C/T(Y)		C/T(Y)
		$callRuleDsRef{"T"}->{"G"}->{"S"}="C";	#	T		G	    C/G(S)		C
		$callRuleDsRef{"T"}->{"G"}->{"M"}="C";	#	T		G	    A/C(M)		C
		$callRuleDsRef{"Y"}->{"G"}->{"Y"}="Y";	#	C/T(Y)	G		C/T(Y)		C/T(Y)
		$callRuleDsRef{"A"}->{"T"}->{"M"}="A";	#	A		T		M			A
		$callRuleDsRef{"A"}->{"T"}->{"W"}="A";	#	A		T		W(A/T)		A
		$callRuleDsRef{"T"}->{"A"}->{"W"}="A";	#	T		A		W(A/T)		T
		$callRuleDsRef{"C"}->{"G"}->{"M"}="C";	#	C		G		M			C
		$callRuleDsRef{"C"}->{"G"}->{"Y"}="C";	#	C		G		C/T(Y)		C
		$callRuleDsRef{"C"}->{"R"}->{"Y"}="Y";	#	C		G/A(R)	C/T(Y)		C/T(Y)
		$callRuleDsRef{"C"}->{"G"}->{"S"}="C";	#	C		G		C/G(S)		C
		$callRuleDsRef{"M"}->{"K"}->{"M"}="M";	#	A/C(M)	K		A/C(M)		A/C(M)
		$callRuleDsRef{"M"}->{"G"}->{"M"}="M";	#	A/C(M)	G		A/C(M)		A/C(M)
		$callRuleDsRef{"M"}->{"T"}->{"M"}="M";	#	A/C(M)	T		A/C(M)		A/C(M)
		$callRuleDsRef{"W"}->{"K"}->{"M"}="M";	#	W		K		A/C(M)		A/C(M)
		$callRuleDsRef{"W"}->{"G"}->{"M"}="M";	#	W		G		A/C(M)		A/C(M)
		$callRuleDsRef{"Y"}->{"K"}->{"M"}="M";	#	C/T(Y)	K		A/C(M)		A/C(M)

my %candidate_snp_info;

sub main(){
	read_snp_list($ARGV[0]);
	load_dbsnp($snp132_file);
	call_genotype();
	report_SNPs();
}
sub report_SNPs(){
	print "SNP position\tSNP call\tSNP qual\tdbSNP\tRefAlleles\tSNP call(fwd)\tAllele count(fwd)\tSNP call(rev)\tAllele count(rev)\n";
	foreach my $index(keys(%candidate_snp_info)){
		next if(!$candidate_snp_info{$index}->{"call"} 
				|| $candidate_snp_info{$index}->{"qual"} < 30
				|| length($candidate_snp_info{$index}->{"call"})>1);
		print "$index\t", $candidate_snp_info{$index}->{"call"}, "\t", $candidate_snp_info{$index}->{"qual"}, "\t";
		if($candidate_snp_info{$index}->{"rs"}){
			print $candidate_snp_info{$index}->{"rs"}, "\t", $candidate_snp_info{$index}->{"allele"}, "\t";
		}else{
			print "-\t-\t";
		}
		if($candidate_snp_info{$index}->{"W"}->{"call"}){
			print $candidate_snp_info{$index}->{"W"}->{"call"}, "\t", $candidate_snp_info{$index}->{"W"}->{"allele_count"}, "\t";
		}else{
			print "-\t-\t";
		}
		if($candidate_snp_info{$index}->{"C"}->{"call"}){
			print $candidate_snp_info{$index}->{"C"}->{"call"}, "\t", $candidate_snp_info{$index}->{"C"}->{"allele_count"}, "\t";
		}else{
			print "-\t-\t";
		}
		print "\n";
	}
}

sub read_snp_list(){
	my $filename = shift;
	open(INFILE, "$filename") || die("Error in opening file $filename\n");
	while(my $line = <INFILE>){
		chop($line);
		my @fields = split(/\t/, $line);
		my $key = $fields[0] . ":" . $fields[1];
		$candidate_snp_info{$key}->{$fields[3]}->{"refBase"} = $fields[3];
		$candidate_snp_info{$key}->{$fields[3]}->{"call"} = $fields[4];
		$candidate_snp_info{$key}->{$fields[3]}->{"qual"} = $fields[5];
		$candidate_snp_info{$key}->{$fields[3]}->{"depth"} = $fields[6];
		$candidate_snp_info{$key}->{"totaldepth"} += $fields[6];
		my $last_pos = scalar(@fields)-1;
		$candidate_snp_info{$key}->{$fields[3]}->{"allele_count"} = join(",", @fields[7..$last_pos]) if(scalar(@fields)>6);
	}
	close(INFILE);
}

sub call_genotype(){
	foreach my $index(keys(%candidate_snp_info)){
		my ($genotype1, $genotype2, $call, $qual);
		next if($candidate_snp_info{$index}->{"totaldepth"}<8);		
		if($candidate_snp_info{$index}->{"W"}->{"call"} && $candidate_snp_info{$index}->{"C"}->{"call"} && $candidate_snp_info{$index}->{"allele"}){
			if($callRuleDsRef{$candidate_snp_info{$index}->{"W"}->{"call"}}->{$candidate_snp_info{$index}->{"C"}->{"call"}}->{$candidate_snp_info{$index}->{"allele"}}){
				$call = $callRuleDsRef{$candidate_snp_info{$index}->{"W"}->{"call"}}->{$candidate_snp_info{$index}->{"C"}->{"call"}}->{$candidate_snp_info{$index}->{"allele"}};				        
			}else{
				$call = "DsRef".$rcTable{$candidate_snp_info{$index}->{"C"}->{"call"}} . ":" .$candidate_snp_info{$index}->{"W"}->{"call"} . "?"
			}			
			$qual = $candidate_snp_info{$index}->{"W"}->{"qual"}+$candidate_snp_info{$index}->{"C"}->{"qual"};			
		}elsif($candidate_snp_info{$index}->{"W"}->{"call"} && $candidate_snp_info{$index}->{"C"}->{"call"}){
			if($callRuleDS{$candidate_snp_info{$index}->{"W"}->{"call"}}->{$candidate_snp_info{$index}->{"C"}->{"call"}}){
				$call = $callRuleDS{$candidate_snp_info{$index}->{"W"}->{"call"}}->{$candidate_snp_info{$index}->{"C"}->{"call"}};
			}else{
				$call = "DS-".$rcTable{$candidate_snp_info{$index}->{"C"}->{"call"}} . ":" . $candidate_snp_info{$index}->{"W"}->{"call"} . "?"
			}			
			$qual = $candidate_snp_info{$index}->{"W"}->{"qual"}+$candidate_snp_info{$index}->{"C"}->{"qual"};			
		}elsif($candidate_snp_info{$index}->{"rs"}){	# for known SNP locations		
			if($candidate_snp_info{$index}->{"W"}->{"call"}){				
				if($callRuleW{$candidate_snp_info{$index}->{"allele"}}->{$candidate_snp_info{$index}->{"W"}->{"call"}}){
					$genotype1 = $callRuleW{$candidate_snp_info{$index}->{"allele"}}->{$candidate_snp_info{$index}->{"W"}->{"call"}};
				}else{
					$genotype1 = "?";
				}				
				$qual = $candidate_snp_info{$index}->{"W"}->{"qual"};
			}					
			if($candidate_snp_info{$index}->{"C"}->{"call"}){
				if($callRuleC{$candidate_snp_info{$index}->{"allele"}}->{$candidate_snp_info{$index}->{"C"}->{"call"}}){
					$genotype2 = $callRuleC{$candidate_snp_info{$index}->{"allele"}}->{$candidate_snp_info{$index}->{"C"}->{"call"}};
				}else{
					$genotype2 = "?";
				}
				$qual += $candidate_snp_info{$index}->{"C"}->{"qual"};				
			}	
			if($genotype1 && $genotype2 && $genotype1 ne $genotype2){
				$call="Ref".$genotype1.":".$genotype2;
			}else{
				$call = $genotype1 ? $genotype1 : $genotype2;
			}
		}
		$candidate_snp_info{$index}->{"call"}=$call if($call);
		$candidate_snp_info{$index}->{"qual"}=$qual if($qual);
	}
}

sub load_dbsnp(){
	my $refSnpFile = shift;
	open(INFILE, "$refSnpFile") || die("Error in opening $refSnpFile!\n");
	while(my $line = <INFILE>){
		my @fields = split(/\t/, $line);
		my $index = $fields[1] . ":" . $fields[3];
		if($candidate_snp_info{$index}){
			next if(!$three2one{$fields[9]} );
			$candidate_snp_info{$index}->{"rs"} = $fields[4];
			$candidate_snp_info{$index}->{"allele"} = $fields[6] eq "+" ? $three2one{$fields[9]} : $rcTable{$three2one{$fields[9]}};
		}
	}
	close(INFILE);
}

main();
