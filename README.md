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



