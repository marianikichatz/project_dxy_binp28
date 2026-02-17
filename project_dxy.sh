#!/bin/bash

# we need to filter the vcf file to get only the SNPs, and remove the indels and multi allelic sites
grep "^#" ProjTaxa.vcf > ProjTaxa.SNP.vcf # keep the header lines
grep -v "^#" ProjTaxa.vcf | awk 'length($4)==1 && length($5)==1 && $5 !~ /,/' >> ProjTaxa.SNP.vcf # keep only the SNPs

# we remove the Naxos2 sample from the vcf file, as we don't need it for the dxy analysis
bcftools view -s ^Naxos2 ProjTaxa.SNP.vcf -o ProjTaxa.SNPonly.vcf

# we zip the vcf file 
bgzip ProjTaxa.SNPonly.vcf

# we index the vcf file
tabix ProjTaxa.SNPonly.vcf.gz

# we set the path to the vcf file for the dxy analysis
VCF=~/Project_binp28/ProjTaxa.SNPonly.vcf.gz

# Statistics per site using vcftools

# quality 
vcftools --gzvcf $VCF --site-quality --out quality_site

# depth
vcftools --gzvcf $VCF --site-mean-depth --out depth_site

# missing data 
vcftools --gzvcf $VCF --missing-site --out missing_site

# maf
vcftools --gzvcf $VCF --freq2 --out max_al  --max-alleles 2


# Statistics per individual using vcftools

# genotype quality
# we use bcftools query to get the GQ for each individual at each site
bcftools query -f '%CHROM\t%POS[\t%GQ]\n' $VCF > quality_ind.lqual

# depth
vcftools --gzvcf $VCF --depth --out depth_ind

# missing data
vcftools --gzvcf $VCF --missing-indv --out missing_ind

# before we decde on the filtering parameters, we can look at the distribution of the quality, depth, missing data and maf 
# to see if there are any outliers or if there are any sites that have a very low quality or depth, or a very high missing data or maf
# we use R to plot the distribution of these statistics and decide on the filtering parameters

# FITERING 

MISS=0.9 # sites with more than 10% of missing data will be removed
QUAL= 30 # sites with a quality score lower than 30 will be removed
MIN_DEPTH=5 # sites with a mean depth lower than 5 will be removed
MAX_DEPTH=27 # sites with a mean depth higher than 27 will be removed
MAF=0.033 # sites with a minor allele frequency lower than 0.033 will be removed

# 1st case -> GQ > 15
GQ=15
vcftools --gzvcf $VCF --remove-indels --max-alleles 2 --maf $MAF --max-missing $MISS --minQ $QUAL --minGQ 15 --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --minDP $MIN_DEPTH --maxDP $--recode --stdout | gzip -c > filtered_gq15.vcf.gz


# 2nd case -> GQ > 10
GQ=10
vcftools --gzvcf $VCF --remove-indels --max-alleles 2 --maf $MAF --max-missing $MISS --minQ $QUAL --minGQ $GQ --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > filtered_gq10.vcf.gz


# we create a populations file for the dxy analysis, with the sample name and the population name separated by a tab
# we have 5 samples for each population, and we have 3 populations: 8N, K and Lesina 
echo -e "8N05240\t8N
8N05890\t8N
8N06612\t8N
8N73248\t8N
8N73604\t8N
K006\tK
K010\tK
K011\tK
K015\tK
K019\tK
Lesina_280\tLesina
Lesina_281\tLesina
Lesina_282\tLesina
Lesina_285\tLesina
Lesina_286\tLesina" > populations.txt

# we need to unzip the filtered vcf file to be able to use it for the dxy analysis, as the dxy script doesn't accept gzipped vcf files
gunzip filtered_gq15.vcf.gz
bgzip filtered_gq15.vcf
# we index the filtered vcf file
tabix filtered_gq15.vcf.gz

# we do the same for the filtered_gq10.vcf.gz file
gunzip filtered_gq10.vcf.gz
bgzip filtered_gq10.vcf
tabix filtered_gq10.vcf.gz

# we run the dxy analysis using pixy, a python script that calculates the dxy between two populations using a vcf file and a populations file
# we will use two different sliding window sizes, 50 kb and 10 kb, to see how the results change with the window size

# for the filtered_gq15.vcf.gz file
pixy --stats dxy --vcf filtered_gq15.vcf.gz --populations populations.txt --window_size 50000 --n_cores 4 --chromosomes 'chr5,chrZ' --output_prefix results_50_15 --bypass_invariant_check --output_folder .

pixy --stats dxy --vcf filtered_gq15.vcf.gz --populations populations.txt --window_size 10000 --n_cores 4 --chromosomes 'chr5,chrZ' --output_prefix results_10_15 --bypass_invariant_check --output_folder .

# for the filtered_gq10.vcf.gz file
pixy --stats dxy --vcf filtered_gq10.vcf.gz --populations populations.txt --window_size 50000 --n_cores 4 --chromosomes 'chr5,chrZ' --output_prefix results_50_10 --bypass_invariant_check --output_folder .

pixy --stats dxy --vcf filtered_gq10.vcf.gz --populations populations.txt --window_size 10000 --n_cores 4 --chromosomes 'chr5,chrZ' --output_prefix results_10_10 --bypass_invariant_check --output_folder .

# we decided to test the sliding window of 100kb for better resolution on the filtered_gq15.vcf.gz file
# as it may give us more informative results than the 50 kb window size
pixy --stats dxy --vcf filtered_gq15.vcf.gz --populations populations.txt --window_size 100000 --n_cores 4 --chromosomes 'chr5,chrZ' --output_prefix results_100_15 --bypass_invariant_check --output_folder .