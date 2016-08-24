#!/usr/bin/perl -w

# Program to count variants per gene and/or pathway from multicalled vcf file
# print files of variant counts for genes and pathways; & per SampleID and file of 'interesting' variants

##IMPROVEMENTS tbc... - chr/pos gene ids?! ##

use strict;
use warnings;

my($path1,$file,$sampleID,$Q,$refbuild,$intgenesfile,$pathsfile,$inhousefile)=@ARGV;

#my $sampleID="4028_hg38";
#my $Q=25;
#my $refbuild="hg38";
#my $path1="/users/nhrg/lustre/SourceBioMito/VCFs";
#my $file=$sampleID."_SourceBio_Q".$Q."_Bcftools13_hg38_noIndel_combined.hg38_multianno.vcf";
#my $intgenesfile="/home/nhrg/WORKING_DATA/InHouse_FastuniqBWAFreebayesPipe_Variants_MAFs_281files.txt.gz";
#my $pathsfile="/home/nhrg/WORKING_DATA/CPDB_pathways_genes.tab";

my $GTfield=0;
my $outfile1 = $path1."/".$sampleID."/".$sampleID."_Q".$Q."_Bcftools13_".$refbuild."_VariantCounts.txt";
my $outfile2 = $path1."/".$sampleID."/".$sampleID."_Q".$Q."_Bcftools13_".$refbuild."_GeneCounts.txt";
my $outfile3 = $path1."/".$sampleID."/".$sampleID."_Q".$Q."_Bcftools13_".$refbuild."_PathwayCounts.txt";
my $outfile4 = $path1."/".$sampleID."/".$sampleID."_Q".$Q."_Bcftools13_".$refbuild."_ProteinAlteringVariants.vcf";
my $outfile5 = $path1."/".$sampleID."/".$sampleID."_Q".$Q."_Bcftools13_".$refbuild."_GenicVariants.vcf";
open(OUT4, ">$outfile4") || die "Cannot open file \"$outfile4\" to write to!\n";
open(OUT5, ">$outfile5") || die "Cannot open file \"$outfile5\" to write to!\n";

## Read pathways into %pathgenes
my %pathgenes=();
my %pathgenecount=();
my %intpaths=();
open INPUT2, $pathsfile or die "Cannot open $pathsfile\n";
my $h=<INPUT2>;
loop2: while (<INPUT2>){
	my $Line2=$_;
	chomp $Line2;
	my @lsplit2=split(/\t/,$Line2);
	my @gsplit=split(/\,/,$lsplit2[3]);
		my $cpath=$lsplit2[0]."_".$lsplit2[1]."_".$lsplit2[2];
	foreach my $g (@gsplit){
		if(!exists $pathgenes{$g}{$cpath}){
			$pathgenes{$g}{$cpath}=0;
			}
		}
	if(!exists $pathgenecount{$cpath}){$pathgenecount{$cpath}=scalar(@gsplit);}	
	if(!exists $intpaths{$cpath}){$intpaths{$cpath}="NO";} #changes to "YES" if path contains gene on the 'interesting' list	
}
close INPUT2;

## Read interesting genes into %intgenes
my %intgenes=();
open INPUT3, $intgenesfile or die "Cannot open $intgenesfile\n";
my $h2=<INPUT3>;
loop3: while (<INPUT3>){
	my $Line2=$_;
	chomp $Line2;
	my @lsplit2=split(/\t/,$Line2);
	if(!exists $intgenes{$lsplit2[4]}){
		$intgenes{$lsplit2[4]}=0;
		}
	if(exists $pathgenes{$lsplit2[4]}){
		foreach my $p (keys %{$pathgenes{$lsplit2[4]}}){$intpaths{$p}="YES";}
		}
}
close INPUT3;

my %inhouse=();
## Read in-house MAF into %inhouse
open INPUT4, $inhousefile or die "Cannot open $inhousefile\n";
my $h4=<INPUT4>;
loop4: while (<INPUT4>){
	my $Line2=$_;
	chomp $Line2;
	my @lsplit2=split(/\t/,$Line2);
	if(!exists $inhouse{$lsplit2[0]}{$lsplit2[1]}{$lsplit2[2]}{$lsplit2[3]}){
		$inhouse{$lsplit2[0]}{$lsplit2[1]}{$lsplit2[2]}{$lsplit2[3]}{'maf'}=$lsplit2[8]; ##or [7] depending on best maf????
		$inhouse{$lsplit2[0]}{$lsplit2[1]}{$lsplit2[2]}{$lsplit2[3]}{'IHCount'}=$lsplit2[4];
		$inhouse{$lsplit2[0]}{$lsplit2[1]}{$lsplit2[2]}{$lsplit2[3]}{'het'}=$lsplit2[5];
		$inhouse{$lsplit2[0]}{$lsplit2[1]}{$lsplit2[2]}{$lsplit2[3]}{'hom'}=$lsplit2[6];
		}
}
close INPUT4;

## Read through vcf and assign counts to %paths and %genes
my %genes=();
my %paths=();
my @ids=();
my %idperCounts=();
open INPUT, $path1."/".$sampleID."/".$file or die "Cannot open $file\n";
loop: while (<INPUT>){
	my $Line=$_;
	chomp $Line;
	## Get sampleIds from header
	if($Line=~/\#CHROM/){
		print OUT4 "$Line\n";
		print OUT5 "$Line\n";
		my @headsplit=split(/\t/,$Line);
		for(my $i=9;$i<scalar(@headsplit);$i++){
			push(@ids,$headsplit[$i]);
			if(!exists $idperCounts{$headsplit[$i]}){
					$idperCounts{$headsplit[$i]}{'Tot'}=0;
					$idperCounts{$headsplit[$i]}{'Het'}=0;
					$idperCounts{$headsplit[$i]}{'Hom'}=0;
					$idperCounts{$headsplit[$i]}{'protAHet'}=0;
					$idperCounts{$headsplit[$i]}{'protAHom'}=0;
					$idperCounts{$headsplit[$i]}{'consvHet'}=0;
					$idperCounts{$headsplit[$i]}{'consvHom'}=0;
					$idperCounts{$headsplit[$i]}{'rareHet'}=0;
					$idperCounts{$headsplit[$i]}{'rareHom'}=0;
					$idperCounts{$headsplit[$i]}{'rareConHet'}=0;
					$idperCounts{$headsplit[$i]}{'rareConHom'}=0;
					$idperCounts{$headsplit[$i]}{'rcdHet'}=0;
					$idperCounts{$headsplit[$i]}{'rcdHom'}=0;
					}
				}#for header sample ids
			foreach my $p (keys %pathgenecount){		
				if(!exists $paths{$p}){
					for(my $i=9;$i<scalar(@headsplit);$i++){
						$paths{$p}{$headsplit[$i]}{'Tot'}=0;
						$paths{$p}{$headsplit[$i]}{'Hom'}=0;
						$paths{$p}{$headsplit[$i]}{'Het'}=0;
						$paths{$p}{$headsplit[$i]}{'protAHet'}=0;
						$paths{$p}{$headsplit[$i]}{'protAHom'}=0;
						$paths{$p}{$headsplit[$i]}{'rcdHet'}=0;
						$paths{$p}{$headsplit[$i]}{'rcdHom'}=0;
						}
					}
				}#foreach path of %pathgenecount
		}#if #CHROM header match
	
	if($Line=~/\#/){next loop;}
	my @lsplit=();
	if($Line=~/^chr/){@lsplit=split(/\t/,$Line);}
	
	## Min var call QUAL filter <= 20
	if($lsplit[5]<=20){next loop;}
	
	## locate the 'GT' info fields, which can differ depending on variant caller
	my @format=split(/\:/,$lsplit[8]);
	for(my $d=0;$d<scalar(@format);$d++){
		if($format[$d]=~/^GT$/){$GTfield=$d;}
		}
	
	## Get gene name and variant type e.g. non-syn, splicing etc.
	my %info=();
	my $gene="NANA"; #any multi-gene name situations??!!!
	my $type="NANA"; #non-syn, syn, splicing, frameshift, etc.
	my $type2="NANA"; #exonic, splicing. intronic, UTR, etc.
	my $rare=0;
	my $vrare=0;
	my $consv=0;
	my $damag=0;
	my $noCD=0;
	my $PA="NO";
	my $GENIC="NO";
	my $rcdhet="NO";
	my $rcdhom="NO";
	my $con="NO";
	my $rarhet="NO";
	my $rarhom="NO";
	my $InH_hets=0;
	my $InH_homs=0;
	
	## Initialise/get values of rare, damaging and conserved measures
	if($lsplit[7]=~/Func.knownGene\=(\S+?);/){$type2=$1;}
	if($lsplit[7]=~/Gene.knownGene\=(\S+?);/){$gene=$1;}
	if($lsplit[7]=~/ExonicFunc.knownGene\=(\S+?);/){$type=$1;}
	if($type=~/^\.$/ or $type=~/^unknown$/){$type=$type2;}
	if($type2=~/^exonic\S+splicing$/){$type="splicing";}
	if($Line=~/\;phyloP20way_mammalian\=(\S+?)\;/){$info{'phyloP20way_mammalian'}=$1;}
	if($Line=~/\;phastCons20way_mammalian\=(\S+?)\;/){$info{'phastCons20way_mammalian'}=$1;}
	if($Line=~/\;SiPhy_29way_logOdds\=(\S+?)\;/){$info{'SiPhy_29way_logOdds'}=$1;}
	if($Line=~/\;SIFT_pred\=(\S+?)\;/){$info{'SIFT_pred'}=$1;}
	if($Line=~/\;Polyphen2_HDIV_pred\=(\S+?)\;/){$info{'Polyphen2_HDIV_pred'}=$1;}
	if($Line=~/\;Polyphen2_HVAR_pred\=(\S+?)\;/){$info{'Polyphen2_HVAR_pred'}=$1;}
	if($Line=~/\;LRT_pred\=(\S+?)\;/){$info{'LRT_pred'}=$1;}
	if($Line=~/\;MutationTaster_pred\=(\S+?)\;/){$info{'MutationTaster_pred'}=$1;}
	if($Line=~/\;MutationAssessor_pred\=(\S+?)\;/){$info{'MutationAssessor_pred'}=$1;}
	if($Line=~/\;FATHMM_pred\=(\S+?)\;/){$info{'FATHMM_pred'}=$1;}
	if($Line=~/\;PROVEAN_pred\=(\S+?)\;/){$info{'PROVEAN_pred'}=$1;}
	if($Line=~/\;VEST3_score\=(\S+?)\;/){$info{'VEST3_score'}=$1;}
	if($Line=~/\;CADD_phred\=(\S+?)\;/){$info{'CADD_phred'}=$1;}
	if($Line=~/\;DANN_score\=(\S+?)\;/){$info{'DANN_score'}=$1;}
	if($Line=~/\;fathmm-MKL_coding_pred\=(\S+?)\;/){$info{'fathmm-MKL_coding_pred'}=$1;}
	if($Line=~/\;MetaSVM_pred\=(\S+?)\;/){$info{'MetaSVM_pred'}=$1;}
	if($Line=~/\;MetaLR_pred\=(\S+?)\;/){$info{'MetaLR_pred'}=$1;}
	if($Line=~/\;integrated_fitCons_score\=(\S+?)\;/){$info{'integrated_fitCons_score'}=$1;}
	if($Line=~/\;ExAC_ALL\=(\S+?)\;/){$info{'ExAC_ALL'}=$1;}
	if($Line=~/\;ExAC_AFR\=(\S+?)\;/){$info{'ExAC_AFR'}=$1;}
	if($Line=~/\;ExAC_FIN\=(\S+?)\;/){$info{'ExAC_FIN'}=$1;}
	if($Line=~/\;esp6500siv2_all\=(\S+?)\;/){$info{'esp6500siv2all'}=$1;}
	if($Line=~/\;HRC_AF\=(\S+?)\;/){$info{'HRC_AF'}=$1;}		
	if($Line=~/\;HRC_non1000G_AF\=(\S+?)\;/){$info{'HRC_non1000G_AF'}=$1;}
	if($Line=~/\;Kaviar_AF\=(\S+?)\;/){$info{'Kaviar_AF'}=$1;}
		
	## Get cumulative rare, conserved and damaging scores
	unless($info{'phyloP20way_mammalian'}=~/^\.$/){if($info{'phyloP20way_mammalian'}>=0.2){$consv++;}}
	unless($info{'phastCons20way_mammalian'}=~/^\.$/){if($info{'phastCons20way_mammalian'}>=0.8){$consv++;}}
	unless($info{'SiPhy_29way_logOdds'}=~/^\.$/){if($info{'SiPhy_29way_logOdds'}>=10){$consv++;}}
	if($info{'phyloP20way_mammalian'}=~/^\.$/ and $info{'phastCons20way_mammalian'}=~/^\.$/ and $info{'SiPhy_29way_logOdds'}=~/^\.$/){$noCD++;}
	if($info{'SIFT_pred'}=~/D/){$damag++;}
	if($info{'Polyphen2_HDIV_pred'}=~/[DP]/){$damag++;}
	if($info{'Polyphen2_HVAR_pred'}=~/[DP]/){$damag++;}
	if($info{'LRT_pred'}=~/D/){$damag++;}
	if($info{'MutationTaster_pred'}=~/D/){$damag++;}
	if($info{'MutationAssessor_pred'}=~/D/){$damag++;}						
	if($info{'FATHMM_pred'}=~/D/){$damag++;}				
	if($info{'PROVEAN_pred'}=~/D/){$damag++;}				
	if($info{'fathmm-MKL_coding_pred'}=~/D/){$damag++;}
	if($info{'MetaSVM_pred'}=~/D/){$damag++;}
	if($info{'MetaLR_pred'}=~/D/){$damag++;}
	unless($info{'VEST3_score'}=~/^\.$/){if($info{'VEST3_score'}>=0.9){$damag++;}}
	unless($info{'CADD_phred'}=~/^\.$/){if($info{'CADD_phred'}>=10){$damag++;}}
	unless($info{'DANN_score'}=~/^\.$/){if($info{'DANN_score'}>=0.9){$damag++;}}
	unless($info{'integrated_fitCons_score'}=~/^\.$/){if($info{'integrated_fitCons_score'}>=0.7){$damag++;}}
	if($info{'SIFT_pred'}=~/\./ and $info{'Polyphen2_HDIV_pred'}=~/\./ and $info{'Polyphen2_HVAR_pred'}=~/\./ and $info{'LRT_pred'}=~/\./ and
		$info{'MutationTaster_pred'}=~/\./ and $info{'MutationAssessor_pred'}=~/\./ and $info{'FATHMM_pred'}=~/\./ and $info{'PROVEAN_pred'}=~/\./ and
			$info{'fathmm-MKL_coding_pred'}=~/\./ and $info{'MetaSVM_pred'}=~/\./ and $info{'MetaLR_pred'}=~/\./ and $info{'VEST3_score'}=~/^\.$/ and 
				$info{'CADD_phred'}=~/^\.$/ and $info{'DANN_score'}=~/^\.$/ and $info{'integrated_fitCons_score'}=~/^\.$/){$noCD++;}
	if($info{'ExAC_ALL'}=~/^\.$/){$rare++;$vrare++;}elsif($info{'ExAC_ALL'}<0.001){$rare++;$vrare++;}elsif($info{'ExAC_ALL'}<0.01){$rare++;}
	if($info{'ExAC_AFR'}=~/^\.$/){$rare++;$vrare++;}elsif($info{'ExAC_AFR'}<0.001){$rare++;$vrare++;}elsif($info{'ExAC_AFR'}<0.01){$rare++;}
	if($info{'ExAC_FIN'}=~/^\.$/){$rare++;$vrare++;}elsif($info{'ExAC_FIN'}<0.001){$rare++;$vrare++;}elsif($info{'ExAC_FIN'}<0.01){$rare++;}
	if($info{'esp6500siv2all'}=~/\.$/){$rare++;$vrare++;}elsif($info{'esp6500siv2all'}<0.001){$rare++;$vrare++;}elsif($info{'esp6500siv2all'}<0.01){$rare++;}
	if($info{'HRC_AF'}=~/^\.$/){$rare++;$vrare++;}elsif($info{'HRC_AF'}<0.001){$rare++;$vrare++;}elsif($info{'HRC_AF'}<0.01){$rare++;}
	if($info{'HRC_non1000G_AF'}=~/^\.$/){$rare++;$vrare++;}elsif($info{'HRC_non1000G_AF'}<0.001){$rare++;$vrare++;}elsif($info{'HRC_non1000G_AF'}<0.01){$rare++;}
	if($info{'Kaviar_AF'}=~/^\.$/){$rare++;$vrare++;}elsif($info{'Kaviar_AF'}<0.001){$rare++;$vrare++;}elsif($info{'Kaviar_AF'}<0.01){$rare++;}
	if(exists $inhouse{$lsplit[0]}{$lsplit[1]}{$lsplit[3]}{$lsplit[4]}){$InH_hets=$inhouse{$lsplit[0]}{$lsplit[1]}{$lsplit[3]}{$lsplit[4]}{'het'};}
	if(exists $inhouse{$lsplit[0]}{$lsplit[1]}{$lsplit[3]}{$lsplit[4]}){$InH_homs=$inhouse{$lsplit[0]}{$lsplit[1]}{$lsplit[3]}{$lsplit[4]}{'hom'};}
		
	## Classify variant type (i.e. protein altering, rare, conserved)
	if($type !~/intergenic/ and $type !~/stream/){$GENIC="YES";} ##!'up/downstream', 'intergenic'
	if($type !~/intergenic/ and $type !~/^synonymous/ and $type !~/intronic/ and $type !~/exonic/ and $type !~/UTR/ and $type !~/stream/){$PA="YES";} ##!'intronic', 'up/downstream', 'exonic', UTR (or splicing) added as label if type is 'unknown' or '.'
	if($consv > 0){$con="YES";}
	if($vrare == 7 and $InH_hets < 5 and $InH_homs < 2){$rarhet="YES";}
	if($rare == 7 and $InH_hets < 5 and $InH_homs < 5){$rarhom="YES";}
	if(($consv > 0 and $vrare == 7 and $damag > 4 and $InH_hets < 5 and $InH_homs < 2 and $type !~/splicing/ and $type !~/[er]tion/) or
		($consv > 0 and $vrare == 7 and $damag > 0 and $InH_hets < 5 and $InH_homs < 2 and $type =~/splicing/) or
		($vrare == 7 and $InH_hets < 5 and $InH_homs < 2 and $type =~/^frameshift_\S+[er]tion/) or
		($vrare == 7 and $InH_hets < 5 and $InH_homs < 2 and $noCD==2 and $type !~/splicing/ and $type !~/[er]tion/ and $type !~/UTR/ and $type !~/stream/)
		){$rcdhet="YES";}
	if(($consv > 0 and $rare == 7 and $damag > 4 and $InH_hets < 5 and $InH_homs < 5 and $type !~/splicing/ and $type !~/[er]tion/) or
		($consv > 0 and $rare == 7 and $damag > 0 and $InH_hets < 5 and $InH_homs < 5 and $type =~/splicing/) or
		($rare == 7 and $InH_hets < 5 and $InH_homs < 5 and $type =~/^frameshift_\S+[er]tion/) or
		($rare == 7 and $InH_hets < 5 and $InH_homs < 5 and $noCD==2 and $type !~/splicing/ and $type !~/[er]tion/ and $type !~/UTR/ and $type !~/stream/)
		){$rcdhom="YES";}
		
	## for multi sample vcfs $ids[$c-9] = match id name to genotype
	for(my $c=9;$c<scalar(@lsplit);$c++){
		my @info=split(/\:/,$lsplit[$c]);
		my $genot=$info[$GTfield];
	
		## Get total counts per sample id of All, het, hom variants
		if($genot!~/0\/0/ and $genot!~/\.\/\./){$idperCounts{$ids[($c-9)]}{'Tot'}++;}
		if($genot=~/0\/[123456789]/ or $genot=~/1\/[23456789]/ or $genot=~/2\/[3456789]/ or $genot=~/3\/[456789]/){$idperCounts{$ids[($c-9)]}{'Het'}++;}
		if($genot=~/1\/1/ or $genot=~/2\/2/ or $genot=~/3\/3/ or $genot=~/4\/4/){$idperCounts{$ids[($c-9)]}{'Hom'}++;}
		
		## Get per sample counts of Exonic & Intronic variants
		if($lsplit[7]=~/Func.knownGene\=exonic/ or $lsplit[7]=~/Func.knownGene\=splicing/ or 
			$lsplit[7]=~/Func.knownGene\=intronic/ or $lsplit[7]=~/Func.knownGene\=\S+?UTR/ or 
			$lsplit[7]=~/Func.knownGene\=\S+?stream/){ #more than one gene !?					
	
		##Initialise per Gene counts ##!!!! will only provide counts for genes present in vcf and not show genes without variation!!!
		if(!exists $genes{$gene}{$ids[($c-9)]}){
			$genes{$gene}{$ids[($c-9)]}{'Tot'}=0;
			$genes{$gene}{$ids[($c-9)]}{'Hom'}=0;
			$genes{$gene}{$ids[($c-9)]}{'Het'}=0;
			$genes{$gene}{$ids[($c-9)]}{'protAHet'}=0;		
			$genes{$gene}{$ids[($c-9)]}{'protAHom'}=0;
			$genes{$gene}{$ids[($c-9)]}{'rcdHet'}=0;		
			$genes{$gene}{$ids[($c-9)]}{'rcdHom'}=0;
			}
		
		#Total Counts for genes and pathways
		if($genot!~/0\/0/ and $genot!~/\.\/\./){
			$genes{$gene}{$ids[($c-9)]}{'Tot'}++;
			if(exists $pathgenes{$gene}){
				foreach my $p (keys %{$pathgenes{$gene}}){
				$paths{$p}{$ids[($c-9)]}{'Tot'}++;
				}}
			}
		
		#Heterozygous Counts	
		if($genot=~/0\/[12345678]/ or $genot=~/1\/[22345678]/ or $genot=~/2\/[345678]/){
			$genes{$gene}{$ids[($c-9)]}{'Het'}++;
			if($PA=~/YES/){
				$idperCounts{$ids[($c-9)]}{'protAHet'}++;
				$genes{$gene}{$ids[($c-9)]}{'protAHet'}++;
				}
			if($con=~/YES/ and $PA=~/YES/){
				$idperCounts{$ids[($c-9)]}{'consvHet'}++;
				}
			if($rarhet=~/YES/ and $PA=~/YES/){
				$idperCounts{$ids[($c-9)]}{'rareHet'}++;
				}
			if($rarhet=~/YES/ and $con=~/YES/ and $PA=~/YES/){
				$idperCounts{$ids[($c-9)]}{'rareConHet'}++;
				}
			if($rcdhet=~/YES/ and $PA=~/YES/){
				$idperCounts{$ids[($c-9)]}{'rcdHet'}++;
				$genes{$gene}{$ids[($c-9)]}{'rcdHet'}++;
				}	
			if(exists $pathgenes{$gene}){
				foreach my $p (keys %{$pathgenes{$gene}}){
					$paths{$p}{$ids[($c-9)]}{'Het'}++;
					if($PA=~/YES/){$paths{$p}{$ids[($c-9)]}{'protAHet'}++;}
					if($rcdhet=~/YES/ and $PA=~/YES/){$paths{$p}{$ids[($c-9)]}{'rcdHet'}++;}
					}}
			}#if heterozygous
		
		#Homozygous Counts
		if($genot=~/1\/1/ or $genot=~/2\/2/ or $genot=~/3\/3/ or $genot=~/4\/4/ or $genot=~/5\/5/){
			$genes{$gene}{$ids[($c-9)]}{'Hom'}++;
			if($PA=~/YES/){
				$idperCounts{$ids[($c-9)]}{'protAHom'}++;
				$genes{$gene}{$ids[($c-9)]}{'protAHom'}++;
				}
			if($con=~/YES/ and $PA=~/YES/){
				$idperCounts{$ids[($c-9)]}{'consvHom'}++;
				}
			if($rarhom=~/YES/ and $PA=~/YES/){
				$idperCounts{$ids[($c-9)]}{'rareHom'}++;
				}
			if($rarhom=~/YES/ and $con=~/YES/ and $PA=~/YES/){
				$idperCounts{$ids[($c-9)]}{'rareConHom'}++;
				}
			if($rcdhom=~/YES/ and $PA=~/YES/){
				$idperCounts{$ids[($c-9)]}{'rcdHom'}++;
				$genes{$gene}{$ids[($c-9)]}{'rcdHom'}++;
				}
			if(exists $pathgenes{$gene}){
				foreach my $p (keys %{$pathgenes{$gene}}){
					$paths{$p}{$ids[($c-9)]}{'Hom'}++;
					if($PA=~/YES/){$paths{$p}{$ids[($c-9)]}{'protAHom'}++;}
					if($rcdhom=~/YES/ and $PA=~/YES/){$paths{$p}{$ids[($c-9)]}{'rcdHom'}++;}
					}}
				}#if homozygous
			}#if within gene exonic/splicing, UTR, up/downstream
		}#for each GT:DP:etc multisample call	
	if($PA =~/YES/){print OUT4 "$Line\n";}
	if($GENIC =~/YES/){print OUT5 "$Line\n";}
}#while INPUT loop
close INPUT;
close OUT4;

## print variant counts to outfiles
open(OUT1, ">$outfile1") || die "Cannot open file \"$outfile1\" to write to!\n";		
print OUT1 "SampleID\tTotal\tHeterozygous\tHomozygous\tProteinAlt\.Het\tProteinAlt\.Hom\tRare\.Het\tRare\.Hom\tConserved\.Het\tConserved\.Hom\tRareConserved\.Het\tRareConserved\.Hom\tRCD\.Het\tRCD\.Hom\n";
foreach my $i (sort keys %idperCounts){		
	print OUT1 "$i\t$idperCounts{$i}{'Tot'}\t$idperCounts{$i}{'Het'}\t$idperCounts{$i}{'Hom'}\t$idperCounts{$i}{'protAHet'}\t$idperCounts{$i}{'protAHom'}\t$idperCounts{$i}{'rareHet'}\t$idperCounts{$i}{'rareHom'}\t$idperCounts{$i}{'consvHet'}\t$idperCounts{$i}{'consvHom'}\t$idperCounts{$i}{'rareConHet'}\t$idperCounts{$i}{'rareConHom'}\t$idperCounts{$i}{'rcdHet'}\t$idperCounts{$i}{'rcdHom'}\n";
}
close OUT1;

open(OUT2, ">$outfile2") || die "Cannot open file \"$outfile2\" to write to!\n";		
print OUT2 "GeneName\tInterest_Gene";
#foreach my $n (@ids){print OUT2 "\t$n\.Total";}
#foreach my $n (@ids){print OUT2 "\t$n\.TotHet";}
#foreach my $n (@ids){print OUT2 "\t$n\.TotHom";}
#foreach my $n (@ids){print OUT2 "\t$n\.ProtAHet";}
#foreach my $n (@ids){print OUT2 "\t$n\.ProtAHom";}
foreach my $n (@ids){print OUT2 "\t$n\.RCDHet";}
foreach my $n (@ids){print OUT2 "\t$n\.RCDHom";}
print OUT2 "\n";
foreach my $g (sort keys %genes){
		print OUT2 "$g";
		if(exists $intgenes{$g}){print OUT2 "\tYES";}
		if(!exists $intgenes{$g}){print OUT2 "\tNO";}
	#foreach my $i (sort keys %{$genes{$g}}){print OUT2 "\t$genes{$g}{$i}{'Tot'}";}
	#foreach my $i (sort keys %{$genes{$g}}){print OUT2 "\t$genes{$g}{$i}{'Het'}";}
	#foreach my $i (sort keys %{$genes{$g}}){print OUT2 "\t$genes{$g}{$i}{'Hom'}";}
	#foreach my $i (sort keys %{$genes{$g}}){print OUT2 "\t$genes{$g}{$i}{'protAHet'}";}
	#foreach my $i (sort keys %{$genes{$g}}){print OUT2 "\t$genes{$g}{$i}{'protAHom'}";}
	foreach my $i (sort keys %{$genes{$g}}){print OUT2 "\t$genes{$g}{$i}{'rcdHet'}";}
	foreach my $i (sort keys %{$genes{$g}}){print OUT2 "\t$genes{$g}{$i}{'rcdHom'}";}
	print OUT2 "\n";
}
close OUT2;

##!!!!Add array lists of sample ids, cases and controls, to separately sum counts and print per category
my @cases=("","");
my @controls=("","");
open(OUT3, ">$outfile3") || die "Cannot open file \"$outfile3\" to write to!\n";
print OUT3 "PathName\tTotal_Genes\tInterest_GenePath";
#foreach my $n (@ids){print OUT3 "\t$n\.Total";}
#foreach my $n (@ids){print OUT3 "\t$n\.TotHet";}
#foreach my $n (@ids){print OUT3 "\t$n\.TotHom";}
#foreach my $n (@ids){print OUT3 "\t$n\.ProtAHet";}
#foreach my $n (@ids){print OUT3 "\t$n\.ProtAHom";}
foreach my $n (@ids){print OUT3 "\t$n\.RCDHet";}
foreach my $n (@ids){print OUT3 "\t$n\.RCDHom";}
print OUT3 "\n";
foreach my $p (sort keys %paths){
		print OUT3 "$p\t$pathgenecount{$p}\t$intpaths{$p}";
	#foreach my $i (sort keys %{$paths{$p}}){print OUT3 "\t$paths{$p}{$i}{'Tot'}";}
	#foreach my $i (sort keys %{$paths{$p}}){print OUT3 "\t$paths{$p}{$i}{'Het'}";}
	#foreach my $i (sort keys %{$paths{$p}}){print OUT3 "\t$paths{$p}{$i}{'Hom'}";}
	#foreach my $i (sort keys %{$paths{$p}}){print OUT3 "\t$paths{$p}{$i}{'protAHet'}";}
	#foreach my $i (sort keys %{$paths{$p}}){print OUT3 "\t$paths{$p}{$i}{'protAHom'}";}
	foreach my $i (sort keys %{$paths{$p}}){print OUT3 "\t$paths{$p}{$i}{'rcdHet'}";}
	foreach my $i (sort keys %{$paths{$p}}){print OUT3 "\t$paths{$p}{$i}{'rcdHom'}";}
	print OUT3 "\n";
}
close OUT3;

exit;
