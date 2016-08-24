#!/usr/bin/perl

#Script to convert vcf file to tab delim text file that can be input to excel

use strict;
use warnings;

my($path1,$file,$sampleID,$intgenesfile,$inhousefile,$omimfile)=@ARGV;

#define variables
my $vcf_file = $path1."/".$sampleID."/".$file;
my $txt_file=$vcf_file."_Total_in_Vars.txt";
my $txt_file_RDC=$vcf_file."_RDC.txt";
my %info;
my %intgenes=();
my %inhouse=();
my %full_genenames=();
my %combined_names=();

## Read interesting genes into %intgenes
open INPUT3, $intgenesfile or die "Cannot open $intgenesfile\n";
my $h3=<INPUT3>;
loop3: while (<INPUT3>){
	my $Line2=$_;
	chomp $Line2;
	my @lsplit2=split(/\t/,$Line2);
	if(!exists $intgenes{$lsplit2[4]}){
		$intgenes{$lsplit2[4]}=0;
		}
}
close INPUT3;

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

#Get full genenames into hash and then combine all terms for same chr-pos/(gene)
open INPUT_F, $omimfile or die "Cannot open $omimfile\n";
	my $head_e =<INPUT_F>;
	loopf: while (<INPUT_F>){
		   my $Line_f=$_;
		   chomp $Line_f;
		   my @linesplit_f = split(/\t/,$Line_f);
		   my $ch=$linesplit_f[0];
		   if($ch!~/^\d/ and $ch!~/^[XY]$/ and $ch!~/^MT$/){next loopf;}##only include numeric chrs ie. not HG183_PATCH
		   if($ch=~/MT/){$ch="M";}
		   my $stpos=$linesplit_f[1];
		   my $endpos=$linesplit_f[2];
		   my $geneN="";
		   if(defined $linesplit_f[3]){$geneN=$linesplit_f[3];}
		   my $ensg="";
		   if(defined $linesplit_f[4]){$ensg=$linesplit_f[4];}
		   my $goT="";
		   if(defined $linesplit_f[5]){$goT=$linesplit_f[5];}
		   my $wikiD="";
		   if(defined $linesplit_f[6]){$wikiD=$linesplit_f[6];}
		   my $mimD="";
		   if(defined $linesplit_f[8]){$mimD=$linesplit_f[8];}
		   if($wikiD=~/\S+/ and !exists $full_genenames{$ch}{$stpos}{$endpos}{'wiki'}{$wikiD}){$full_genenames{$ch}{$stpos}{$endpos}{'wiki'}{$wikiD}=0;}
		   if($goT=~/\S+/ and !exists $full_genenames{$ch}{$stpos}{$endpos}{'go'}{$goT}){$full_genenames{$ch}{$stpos}{$endpos}{'go'}{$goT}=0;}
		   if($mimD=~/\S+/ and !exists $full_genenames{$ch}{$stpos}{$endpos}{'mim'}{$mimD}){$full_genenames{$ch}{$stpos}{$endpos}{'mim'}{$mimD}=0;}
		   if($geneN=~/\S+/ and !exists $full_genenames{$ch}{$stpos}{$endpos}{'gene'}{$geneN}){$full_genenames{$ch}{$stpos}{$endpos}{'gene'}{$geneN}=0;}
		   if($ensg=~/\S+/ and !exists $full_genenames{$ch}{$stpos}{$endpos}{'ensg'}{$ensg}){$full_genenames{$ch}{$stpos}{$endpos}{'ensg'}{$ensg}=0;}
}
close INPUT_F;
foreach my $c (keys %full_genenames){
	foreach my $sp (keys %{$full_genenames{$c}}){
		foreach my $ep (keys %{$full_genenames{$c}{$sp}}){
				
		if(exists $full_genenames{$c}{$sp}{$ep}{'gene'}){	
	foreach my $n (keys %{$full_genenames{$c}{$sp}{$ep}{'gene'}}){
		if(exists $combined_names{$c}{$sp}{$ep}{'gene'}){$combined_names{$c}{$sp}{$ep}{'gene'}=$combined_names{$c}{$sp}{$ep}{'gene'}.";".$n;}
		if(!exists $combined_names{$c}{$sp}{$ep}{'gene'}){$combined_names{$c}{$sp}{$ep}{'gene'}=$n;}
		}}
		if(exists $full_genenames{$c}{$sp}{$ep}{'ensg'}){	
	foreach my $n (keys %{$full_genenames{$c}{$sp}{$ep}{'ensg'}}){
		if(exists $combined_names{$c}{$sp}{$ep}{'ensg'}){$combined_names{$c}{$sp}{$ep}{'ensg'}=$combined_names{$c}{$sp}{$ep}{'ensg'}.";".$n;}
		if(!exists $combined_names{$c}{$sp}{$ep}{'ensg'}){$combined_names{$c}{$sp}{$ep}{'ensg'}=$n;}
		}}
		if(exists $full_genenames{$c}{$sp}{$ep}{'wiki'}){	
	foreach my $w (keys %{$full_genenames{$c}{$sp}{$ep}{'wiki'}}){
		if(exists $combined_names{$c}{$sp}{$ep}{'wiki'}){$combined_names{$c}{$sp}{$ep}{'wiki'}=$combined_names{$c}{$sp}{$ep}{'wiki'}.";".$w;}
		if(!exists $combined_names{$c}{$sp}{$ep}{'wiki'}){$combined_names{$c}{$sp}{$ep}{'wiki'}=$w;}
		}}
		if(exists $full_genenames{$c}{$sp}{$ep}{'go'}){
	foreach my $g (keys %{$full_genenames{$c}{$sp}{$ep}{'go'}}){
		if(exists $combined_names{$c}{$sp}{$ep}{'go'}){$combined_names{$c}{$sp}{$ep}{'go'}=$combined_names{$c}{$sp}{$ep}{'go'}.";".$g;}
		if(!exists $combined_names{$c}{$sp}{$ep}{'go'}){$combined_names{$c}{$sp}{$ep}{'go'}=$g;}
		}}
		if(exists $full_genenames{$c}{$sp}{$ep}{'mim'}){
	foreach my $m (keys %{$full_genenames{$c}{$sp}{$ep}{'mim'}}){
		if(exists $combined_names{$c}{$sp}{$ep}{'mim'}){$combined_names{$c}{$sp}{$ep}{'mim'}=$combined_names{$c}{$sp}{$ep}{'mim'}.";".$m;}
		if(!exists $combined_names{$c}{$sp}{$ep}{'mim'}){$combined_names{$c}{$sp}{$ep}{'mim'}=$m;}
		}}
}}}

#Read each VCF file line into %info and print to txt files
open(OUT, ">$txt_file") || die "Cannot open file \"$txt_file\" to write to!\n";
open(OUT2, ">$txt_file_RDC") || die "Cannot open file \"$txt_file_RDC\" to write to!\n";
open INPUT2, $vcf_file or die "Cannot open $vcf_file\n";
	loop2: while (<INPUT2>){
		my $Line=$_;
		chomp $Line;
		if($Line=~/\#\#/){next loop2;}
		my @lsplit=split(/\t/,$Line);
		
		#extract header line info and print relevant column headers to out file
		if($Line=~/\#CHROM/){
			#SAMBCF ##CHROM	POS	ID	REF	ALT	QUAL	FILTER	(INFO1)	FORMAT	NDCN_1175_AGGTAAGA_L005	NDCN_2453_TCTAGCGT_L006	WTCCC125661_AGTAGATC_L003	WTCCC125662_AGACTATA_L003	WTCCC126287_TCACAGCA_L003	WTCCC126463_CACGAGAT_L003
			for(my $c=0;$c<7;$c++){print OUT "$lsplit[$c]\t"; print OUT2 "$lsplit[$c]\t";}
			my @headA=("DP","DP4","AN","AC","MQ","InterestGene",
			"FuncknownGene","GeneknownGene","GeneDetailknownGene","ExonicFuncknownGene","AAChangeknownGene",
			"InhouseExN","InhouseHets","InhouseHoms","Rare0.01_Count7","vRare0.001_Count7","ExAC_ALL","ExAC_AFR","ExAC_FIN",
			"esp6500siv2all","HRC_AF","HRC_non1000G_AF","Kaviar_AF","avsnp144",
			"clinvar20151201","Damage_Count15","SIFT_pred","Polyphen2_HDIV_pred","Polyphen2_HVAR_pred","LRT_pred","MutationTaster_pred","MutationAssessor_pred",
			"FATHMM_pred","PROVEAN_pred","VEST3_score","CADD_phred","DANN_score","fathmm-MKL_coding_pred","MetaSVM_pred","MetaLR_pred",
			"integrated_fitCons_score","GERP++_RS","phyloP20way_mammalian","phastCons20way_mammalian","SiPhy_29way_logOdds","Conservation_Count3",
			"Interpro_domain","dbscSNV_ADA_SCORE","dbscSNV_RF_SCORE","ShortGeneNames","ENSGGeneNames","FullGeneNames","GO-terms","OMIM");
			for(my $c1=0;$c1<scalar(@headA);$c1++){print OUT "$headA[$c1]\t"; print OUT2 "$headA[$c1]\t";}
			for(my $c2=8;$c2<scalar(@lsplit);$c2++){print OUT "$lsplit[$c2]\t"; print OUT2 "$lsplit[$c2]\t";} #genotype only header
			#for(my $c2=8;$c2<scalar(@lsplit);$c2++){print OUT "$lsplit[$c2]\t";}#full geno + info header
			print OUT "END\n";
			print OUT2 "END\n";
			next loop2;
			}#end header line
		
#Initialise %info keys
$info{'DP'}=".";
$info{'DP4'}=".";
$info{'AN'}=".";
$info{'AC'}=".";
$info{'MQ'}=".";
$info{'IntGene'}=".";
$info{'FuncknownGene'}=".";
$info{'GeneknownGene'}=".";
$info{'GeneDetailknownGene'}=".";
$info{'ExonicFuncknownGene'}=".";
$info{'AAChangeknownGene'}=".";
$info{'InH_ExN'}=0;
$info{'InH_hets'}=0;
$info{'InH_homs'}=0;
$info{'InH_MAF'}=0;
$info{'Rare_C'}=0;
$info{'vRare_C'}=0;
$info{'ExAC_ALL'}=".";
$info{'ExAC_AFR'}=".";
$info{'ExAC_FIN'}=".";
$info{'esp6500siv2all'}=".";
$info{'avsnp144'}=".";
$info{'clinvar20151201'}=".";
$info{'Damag_C'}=0;
$info{'SIFT_pred'}=".";
$info{'Polyphen2_HDIV_pred'}=".";
$info{'Polyphen2_HVAR_pred'}=".";
$info{'LRT_pred'}=".";
$info{'MutationTaster_pred'}=".";
$info{'MutationAssessor_pred'}=".";
$info{'FATHMM_pred'}=".";
$info{'PROVEAN_pred'}=".";
$info{'VEST3_score'}=".";
$info{'CADD_phred'}=".";
$info{'DANN_score'}=".";
$info{'fathmm-MKL_coding_pred'}=".";
$info{'MetaSVM_pred'}=".";
$info{'MetaLR_pred'}=".";
$info{'integrated_fitCons_score'}=".";
$info{'integrated_confidence_value'}=".";
$info{'GERP++_RS'}=".";
$info{'Conserve_C'}=0;
$info{'phyloP20way_mammalian'}=".";
$info{'phastCons20way_mammalian'}=".";
$info{'SiPhy_29way_logOdds'}=".";
$info{'Interpro_domain'}=".";
$info{'dbscSNV_ADA_SCORE'}=".";
$info{'dbscSNV_RF_SCORE'}=".";
$info{'HRC_AF'}=".";
$info{'HRC_non1000G_AF'}=".";
$info{'Kaviar_AF'}=".";
$info{'nci60'}=".";
$info{'ShortGeneNames'}=".";
$info{'ENSGGeneNames'}=".";
$info{'FullGeneNames'}=".";
$info{'GO-terms'}=".";
$info{'OMIM'}=".";

		#extract annotation information
		if($Line=~/DP\=(\S+?)\;/){$info{'DP'}=$1;}
		if($Line=~/\;DP4\=(\S+?)\;/){$info{'DP4'}=$1;}
		if($Line=~/\;AN\=(\S+?)\;/){$info{'AN'}=$1;}
		if($Line=~/\;AC\=(\S+?)\;/){$info{'AC'}=$1;}
		if($Line=~/\;MQ\=(\S+?)\;/){$info{'MQ'}=$1;}
		if($Line=~/\;Func.knownGene\=(\S+?)\;/){$info{'FuncknownGene'}=$1;}
		if($Line=~/\;Gene.knownGene\=(\S+?)\;/){$info{'GeneknownGene'}=$1;}
		if(exists $intgenes{$info{'GeneknownGene'}}){$info{'IntGene'}="YES";}
		if($Line=~/\;GeneDetail.knownGene\=(\S+?)\;/){$info{'GeneDetailknownGene'}=$1;}
		if($Line=~/\;ExonicFunc.knownGene\=(\S+?)\;/){$info{'ExonicFuncknownGene'}=$1;}
		if($Line=~/\;AAChange.knownGene\=(\S+?)\;/){$info{'AAChangeknownGene'}=$1;}
		if(exists $inhouse{$lsplit[0]}{$lsplit[1]}{$lsplit[3]}{$lsplit[4]}){$info{'InH_ExN'}=$inhouse{$lsplit[0]}{$lsplit[1]}{$lsplit[3]}{$lsplit[4]}{'IHCount'};}
		if(exists $inhouse{$lsplit[0]}{$lsplit[1]}{$lsplit[3]}{$lsplit[4]}){$info{'InH_hets'}=$inhouse{$lsplit[0]}{$lsplit[1]}{$lsplit[3]}{$lsplit[4]}{'het'};}
		if(exists $inhouse{$lsplit[0]}{$lsplit[1]}{$lsplit[3]}{$lsplit[4]}){$info{'InH_homs'}=$inhouse{$lsplit[0]}{$lsplit[1]}{$lsplit[3]}{$lsplit[4]}{'hom'};}
		if(exists $inhouse{$lsplit[0]}{$lsplit[1]}{$lsplit[3]}{$lsplit[4]}){$info{'InH_MAF'}=$inhouse{$lsplit[0]}{$lsplit[1]}{$lsplit[3]}{$lsplit[4]}{'maf'};}
		if($Line=~/\;ExAC_ALL\=(\S+?)\;/){$info{'ExAC_ALL'}=$1;}
		if($Line=~/\;ExAC_AFR\=(\S+?)\;/){$info{'ExAC_AFR'}=$1;}
		if($Line=~/\;ExAC_FIN\=(\S+?)\;/){$info{'ExAC_FIN'}=$1;}
		if($Line=~/\;esp6500siv2_all\=(\S+?)\;/){$info{'esp6500siv2all'}=$1;}
		if($Line=~/\;avsnp144\=(\S+?)\;/){$info{'avsnp144'}=$1;}
		if($Line=~/\;clinvar20151201\=(\S+?)\;/){$info{'clinvar20151201'}=$1;}
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
		if($Line=~/\;GERP++_RS\=(\S+?)\;/){$info{'GERP++_RS'}=$1;}
		if($Line=~/\;phyloP20way_mammalian\=(\S+?)\;/){$info{'phyloP20way_mammalian'}=$1;}
		if($Line=~/\;phastCons20way_mammalian\=(\S+?)\;/){$info{'phastCons20way_mammalian'}=$1;}
		if($Line=~/\;SiPhy_29way_logOdds\=(\S+?)\;/){$info{'SiPhy_29way_logOdds'}=$1;}
		if($Line=~/\;Interpro_domain\=(\S+?)\;/){$info{'Interpro_domain'}=$1;}
		if($Line=~/\;dbscSNV_ADA_SCORE\=(\S+?)\;/){$info{'dbscSNV_ADA_SCORE'}=$1;}
		if($Line=~/\;dbscSNV_RF_SCORE\=(\S+?)\;/){$info{'dbscSNV_RF_SCORE'}=$1;}
		if($Line=~/\;HRC_AF\=(\S+?)\;/){$info{'HRC_AF'}=$1;}		
		if($Line=~/\;HRC_non1000G_AF\=(\S+?)\;/){$info{'HRC_non1000G_AF'}=$1;}
		if($Line=~/\;Kaviar_AF\=(\S+?)\;/){$info{'Kaviar_AF'}=$1;}
		if($Line=~/\;nci60\=(\S+?)\;/){$info{'nci60'}=$1;}

		##Get cumulative damage count
		if($info{'SIFT_pred'}=~/D/){$info{'Damag_C'}++;}
		if($info{'Polyphen2_HDIV_pred'}=~/[DP]/){$info{'Damag_C'}++;}
		if($info{'Polyphen2_HVAR_pred'}=~/[DP]/){$info{'Damag_C'}++;}
		if($info{'LRT_pred'}=~/D/){$info{'Damag_C'}++;}
		if($info{'MutationTaster_pred'}=~/D/){$info{'Damag_C'}++;}
		if($info{'MutationAssessor_pred'}=~/D/){$info{'Damag_C'}++;}						
		if($info{'FATHMM_pred'}=~/D/){$info{'Damag_C'}++;}				
		if($info{'PROVEAN_pred'}=~/D/){$info{'Damag_C'}++;}				
		if($info{'fathmm-MKL_coding_pred'}=~/D/){$info{'Damag_C'}++;}
		if($info{'MetaSVM_pred'}=~/D/){$info{'Damag_C'}++;}
		if($info{'MetaLR_pred'}=~/D/){$info{'Damag_C'}++;}
		unless($info{'VEST3_score'}=~/^\.$/){if($info{'VEST3_score'}>=0.9){$info{'Damag_C'}++;}}
		unless($info{'CADD_phred'}=~/^\.$/){if($info{'CADD_phred'}>=10){$info{'Damag_C'}++;}}
		unless($info{'DANN_score'}=~/^\.$/){if($info{'DANN_score'}>=0.9){$info{'Damag_C'}++;}}
		unless($info{'integrated_fitCons_score'}=~/^\.$/){if($info{'integrated_fitCons_score'}>=0.7){$info{'Damag_C'}++;}}		
		
		##Get cumulative conserved amino acid count
		unless($info{'phyloP20way_mammalian'}=~/^\.$/){if($info{'phyloP20way_mammalian'}>=0.2){$info{'Conserve_C'}++;}}
		unless($info{'phastCons20way_mammalian'}=~/^\.$/){if($info{'phastCons20way_mammalian'}>=0.8){$info{'Conserve_C'}++;}}
		unless($info{'SiPhy_29way_logOdds'}=~/^\.$/){if($info{'SiPhy_29way_logOdds'}>=10){$info{'Conserve_C'}++;}}

		##Get cumulative rare count MAF<0.01
		if($info{'ExAC_ALL'}=~/^\.$/){$info{'Rare_C'}++;}elsif($info{'ExAC_ALL'}<0.01){$info{'Rare_C'}++;} 
		if($info{'ExAC_AFR'}=~/^\.$/){$info{'Rare_C'}++;}elsif($info{'ExAC_AFR'}<0.01){$info{'Rare_C'}++;}
		if($info{'ExAC_FIN'}=~/^\.$/){$info{'Rare_C'}++;}elsif($info{'ExAC_FIN'}<0.01){$info{'Rare_C'}++;}
		if($info{'esp6500siv2all'}=~/\.$/){$info{'Rare_C'}++;}elsif($info{'esp6500siv2all'}<0.01){$info{'Rare_C'}++;}
		if($info{'HRC_AF'}=~/^\.$/){$info{'Rare_C'}++;}elsif($info{'HRC_AF'}<0.01){$info{'Rare_C'}++;}
		if($info{'HRC_non1000G_AF'}=~/^\.$/){$info{'Rare_C'}++;}elsif($info{'HRC_non1000G_AF'}<0.01){$info{'Rare_C'}++;}
		if($info{'Kaviar_AF'}=~/^\.$/){$info{'Rare_C'}++;}elsif($info{'Kaviar_AF'}<0.01){$info{'Rare_C'}++;}
		
		##Get cumulative very rare count MAF<0.001
		if($info{'ExAC_ALL'}=~/^\.$/){$info{'vRare_C'}++;}elsif($info{'ExAC_ALL'}<0.001){$info{'vRare_C'}++;} 
		if($info{'ExAC_AFR'}=~/^\.$/){$info{'vRare_C'}++;}elsif($info{'ExAC_AFR'}<0.001){$info{'vRare_C'}++;}
		if($info{'ExAC_FIN'}=~/^\.$/){$info{'vRare_C'}++;}elsif($info{'ExAC_FIN'}<0.001){$info{'vRare_C'}++;}
		if($info{'esp6500siv2all'}=~/^\.$/){$info{'vRare_C'}++;}elsif($info{'esp6500siv2all'}<0.001){$info{'vRare_C'}++;}
		if($info{'HRC_AF'}=~/^\.$/){$info{'vRare_C'}++;}elsif($info{'HRC_AF'}<0.001){$info{'vRare_C'}++;}
		if($info{'HRC_non1000G_AF'}=~/\.$/){$info{'vRare_C'}++;}elsif($info{'HRC_non1000G_AF'}<0.001){$info{'vRare_C'}++;}
		if($info{'Kaviar_AF'}=~/^\.$/){$info{'vRare_C'}++;}elsif($info{'Kaviar_AF'}<0.001){$info{'vRare_C'}++;}
		
		##Get Other GO, phenotype, GeneName annotation
		my $chr=$lsplit[0];
		$chr=~s/chr//;
		my $pos_st=$lsplit[1];
		foreach my $st (sort {$a<=>$b} keys %{$combined_names{$chr}}){
			foreach my $ed (sort {$a<=>$b} keys %{$combined_names{$chr}{$st}}){
				if($pos_st>=$st and $pos_st<=$ed){
					if(exists $combined_names{$chr}{$st}{$ed}{'wiki'}){$info{'FullGeneNames'}=$info{'FullGeneNames'}.":".$combined_names{$chr}{$st}{$ed}{'wiki'};}
					if(exists $combined_names{$chr}{$st}{$ed}{'go'}){$info{'GO-terms'}=$info{'GO-terms'}.":".$combined_names{$chr}{$st}{$ed}{'go'};}
					if(exists $combined_names{$chr}{$st}{$ed}{'mim'}){$info{'OMIM'}=$info{'OMIM'}.":".$combined_names{$chr}{$st}{$ed}{'mim'};}
					if(exists $combined_names{$chr}{$st}{$ed}{'ensg'}){$info{'ENSGGeneNames'}=$info{'ENSGGeneNames'}.":".$combined_names{$chr}{$st}{$ed}{'ensg'};}
					if(exists $combined_names{$chr}{$st}{$ed}{'gene'}){$info{'ShortGeneNames'}=$info{'ShortGeneNames'}.":".$combined_names{$chr}{$st}{$ed}{'gene'};}					
					}
			}}
		
		#print line info to output txt file
		for(my $c=0;$c<7;$c++){print OUT "$lsplit[$c]\t";}
		my @lineA=("$info{'DP'}","$info{'DP4'}","$info{'AN'}","$info{'AC'}","$info{'MQ'}","$info{'IntGene'}",
		"$info{'FuncknownGene'}","$info{'GeneknownGene'}","$info{'GeneDetailknownGene'}","$info{'ExonicFuncknownGene'}","$info{'AAChangeknownGene'}",
		"$info{'InH_ExN'}","$info{'InH_hets'}","$info{'InH_homs'}",
		"$info{'Rare_C'}","$info{'vRare_C'}","$info{'ExAC_ALL'}","$info{'ExAC_AFR'}","$info{'ExAC_FIN'}","$info{'esp6500siv2all'}",
		"$info{'HRC_AF'}","$info{'HRC_non1000G_AF'}","$info{'Kaviar_AF'}",
		"$info{'avsnp144'}","$info{'clinvar20151201'}","$info{'Damag_C'}","$info{'SIFT_pred'}","$info{'Polyphen2_HDIV_pred'}","$info{'Polyphen2_HVAR_pred'}",
		"$info{'LRT_pred'}","$info{'MutationTaster_pred'}","$info{'MutationAssessor_pred'}","$info{'FATHMM_pred'}","$info{'PROVEAN_pred'}",
		"$info{'VEST3_score'}","$info{'CADD_phred'}","$info{'DANN_score'}","$info{'fathmm-MKL_coding_pred'}","$info{'MetaSVM_pred'}",
		"$info{'MetaLR_pred'}","$info{'integrated_fitCons_score'}","$info{'GERP++_RS'}","$info{'phyloP20way_mammalian'}",
		"$info{'phastCons20way_mammalian'}","$info{'SiPhy_29way_logOdds'}","$info{'Conserve_C'}","$info{'Interpro_domain'}","$info{'dbscSNV_ADA_SCORE'}","$info{'dbscSNV_RF_SCORE'}",
		"$info{'ShortGeneNames'}","$info{'ENSGGeneNames'}","$info{'FullGeneNames'}","$info{'GO-terms'}","$info{'OMIM'}");
		for(my $c1=0;$c1<scalar(@lineA);$c1++){print OUT "$lineA[$c1]\t";}
		for(my $c2=8;$c2<scalar(@lsplit);$c2++){my @genosplit=split(/:/,$lsplit[$c2]); print OUT "\'$genosplit[0]\t";}
		#for(my $c2=8;$c2<scalar(@lsplit);$c2++){print OUT "$lsplit[$c2]\t";}
		print OUT "END\n";
		
		# If Rare, conserved and damaging print line info to outfile2
		if(($info{'Rare_C'}==7 and $info{'Conserve_C'}>0 and $info{'Damag_C'}>4 and $info{'InH_hets'}<5 and $info{'InH_homs'}<5) or
			($info{'Rare_C'}==7 and $info{'Conserve_C'}>0 and $info{'Damag_C'}>0 and $info{'InH_hets'}<5 and $info{'InH_homs'}<5 and $info{'FuncknownGene'}=~/splicing/) or
			($info{'Rare_C'}==7 and $info{'InH_hets'}<5 and $info{'InH_homs'}<5 and $info{'ExonicFuncknownGene'}=~/[dn][es][le][er]tion/) or
			($info{'Rare_C'}==7 and $info{'Conserve_C'}==0 and $info{'Damag_C'}==0 and $info{'InH_hets'}<5 and $info{'InH_homs'}<5)
			){
		for(my $c=0;$c<7;$c++){print OUT2 "$lsplit[$c]\t";}
		my @lineA=("$info{'DP'}","$info{'DP4'}","$info{'AN'}","$info{'AC'}","$info{'MQ'}","$info{'IntGene'}",
		"$info{'FuncknownGene'}","$info{'GeneknownGene'}","$info{'GeneDetailknownGene'}","$info{'ExonicFuncknownGene'}","$info{'AAChangeknownGene'}",
		"$info{'InH_ExN'}","$info{'InH_hets'}","$info{'InH_homs'}",
		"$info{'Rare_C'}","$info{'vRare_C'}","$info{'ExAC_ALL'}","$info{'ExAC_AFR'}","$info{'ExAC_FIN'}","$info{'esp6500siv2all'}",
		"$info{'HRC_AF'}","$info{'HRC_non1000G_AF'}","$info{'Kaviar_AF'}",
		"$info{'avsnp144'}","$info{'clinvar20151201'}","$info{'Damag_C'}","$info{'SIFT_pred'}","$info{'Polyphen2_HDIV_pred'}","$info{'Polyphen2_HVAR_pred'}",
		"$info{'LRT_pred'}","$info{'MutationTaster_pred'}","$info{'MutationAssessor_pred'}","$info{'FATHMM_pred'}","$info{'PROVEAN_pred'}",
		"$info{'VEST3_score'}","$info{'CADD_phred'}","$info{'DANN_score'}","$info{'fathmm-MKL_coding_pred'}","$info{'MetaSVM_pred'}",
		"$info{'MetaLR_pred'}","$info{'integrated_fitCons_score'}","$info{'GERP++_RS'}","$info{'phyloP20way_mammalian'}",
		"$info{'phastCons20way_mammalian'}","$info{'SiPhy_29way_logOdds'}","$info{'Conserve_C'}","$info{'Interpro_domain'}","$info{'dbscSNV_ADA_SCORE'}","$info{'dbscSNV_RF_SCORE'}",
		"$info{'ShortGeneNames'}","$info{'ENSGGeneNames'}","$info{'FullGeneNames'}","$info{'GO-terms'}","$info{'OMIM'}");
		for(my $c1=0;$c1<scalar(@lineA);$c1++){print OUT2 "$lineA[$c1]\t";}
		for(my $c2=8;$c2<scalar(@lsplit);$c2++){my @genosplit=split(/:/,$lsplit[$c2]); print OUT2 "\'$genosplit[0]\t";}
		#for(my $c2=8;$c2<scalar(@lsplit);$c2++){print OUT2 "$lsplit[$c2]\t";}
		print OUT2 "END\n";
		}# if RCD
}
close INPUT2;
close OUT;
close OUT2;
exit;
