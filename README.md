# WGS_Intra-locus-sexual-conflict
Script allowing to replicates analysis from the manuscript "Searching for intra-locus sexual conflicts by the mean of whole genome sequencing in the three-spined stickleback (Gasterosteus acculeatus)"


# Setup:
in 02_maps/sex/ you should have a tab-separated file named sexmap
format: ind1\tsex\n ...

as well as a M and F file respectively containing list of male and female samples
in  03_genome
you should have your reference genome (including Y if pertinent)

in 04_datasets/Full_genome you should put your filtered vcf
Recomended format: vcf.gz + bcf (for lostruct)


# Detection of large structural variant
If this manuscript, first PCA analysis rised suspition of potential inversions segregating in our population
To identidy and isolate those region, we relied on lostruc packages, that perform local PCA analysis accross the genome,
Using the 01_Explore_inversions.R 
A lot of parameters (especially plots and sample size) are specific to our study case so pay attention if you want to use it

After this, you can use bcftools to remove identified regions for exemple:
bcftools view -t ^9:3850000-8600000,16:17100000-18050000,21:1800000-2700000 -Oz -o 04_datasets/Full_genome/NO_INV_COV_MISS_QUAL_MAF010_filtered_99_autosomes_biallelic_SOR_sorted_renamed_NoUn.vcf.gz 04_datasets/Full_genome/COV_MISS_QUAL_MAF010_filtered_99_autosomes_biallelic_SOR_sorted_renamed_NoUn.bcf

and use that new VCF for the rest of the pipeline

# Population structure:
