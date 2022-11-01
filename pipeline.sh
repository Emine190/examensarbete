
#First create a list of all the samples or data and have it with the pointer sample so it can be called with $sample command for the loop.
# this sample was 
#Fix the sample.txt file to include the new control women so the program can run and do bioinformatic things. 

samples_to_include=/media/god/data1/Emil/Control/sample

cat $samples_to_include | while read sample
do

echo $sample
control=/media/god/LaCie/GTEx_-50_females_whole_genome_seq/batch2
Input1=/media/god/data1/Emil/Control/untrimmed/$sample'_'1.fastq.gz
Input2=/media/god/data1/Emil/Control/untrimmed/$sample'_'2.fastq.gz

#Have samtools view function where the genome is and becomes a sorted bamfile.
samtools view -b -@32 $control/$sample.cram | samtools sort -o $sample.bam -@32 -n - 

samtools fastq -@32 -1 $Input1 -2 $Input2  $sample.bam

rm $sample.bam

#removes the bam file intermediates so they dont take up unneccesary space on the disks.
# Here pointers for the reference genome and the pointers for the BWA and fastq.
ref_bwa_mem=/media/god/data1/genome_indexes/bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
trimmedin1=/media/god/data1/Emil/Control/trimmed/trimmed.$sample'_'1.fastq.gz
trimmedin2=/media/god/data1/Emil/Control/trimmed/trimmed.$sample'_'2.fastq.gz
ID=$(echo $sample | cut -d'-' -f1,2)
SM=$(echo $sample | cut -d'-' -f1,2)
Aligned=/media/god/data1/Emil/Control/aligned/$sample.aligned.bam
echo $ID
echo $SM
echo $Input1
echo $Input2

# software
sambamba=/media/god/data1/software/sambamba-0.7.1-linux-static
fastp=/media/god/data1/software/fastp
gatk=/media/god/data1/software/gatk-4.2.6.1/gatk
igv=/media/god/data1/software/IGV_Linux_2.8.0

# quality trimming
$fastp -i $Input1 -I $Input2 -o $trimmedin1 -O $trimmedin2

# trimming of the data to remove primers from the sequences.
rm $Input1
rm $Input2
#burrows-Wheeler Alignment tool to align the trimmed data from fastp. 
#the Functions will not work if you take out the X chromosome before the BWA it needs to be done after.

#Program: bwa (alignment via Burrows-Wheeler transformation)
#Version: 0.7.17-r1188

bwa mem -t32 -M -R "@RG\tID:$ID\tSM:$SM\tPL:ILLUMINA:150" $ref_bwa_mem $trimmedin1 $trimmedin2 | samtools sort -o $Aligned -@32 -

#Run an index to run the view functions to create just the X chromosome file.
#Creates the pointer for the index file for the index function so the bai files get to the bam files.
rm $trimmedin1
rm $trimmedin2

bai=/media/god/data1/Emil/Control/aligned/$sample.aligned.bam.bai
#Rerun all the control samples from here in the control group. Something didnt go correctly from here. The pointer to bai was wrong.
samtools index -b -@32 $Aligned $bai

AlignedX=/media/god/data1/Emil/Control/alignedXchr/$sample.AlignedXChr.bam

#This is the point where the aligned genome can be split to just X chromosome
samtools view -b -@32 $Aligned chrX > $AlignedX

#Reverse all the UPIC-tags to the sample tag so the script can flow.
NodupX=/media/god/data1/Emil/Control/NodupXchr/$sample.NodupXchr.bam 
# Finish this pointer with the folder. 
#Samabamba version is 0.7.1-linux-static 
# -r to remove the duplicates located on the file. -r is the function that removes the dups instead of just marking them. 
$sambamba markdup -r -t32 $AlignedX $NodupX
Haplotype=/media/god/data1/Emil/Control/haplotypecaller/$sample.Haplotype.g.vcf.gz
rm $Aligned
#Got a separate version from the UPIC from björn just have to do the indexing.
baiNoX=/media/god/data1/Emil/Control/NodupXchr/$sample.NodupXchr.bam.bai
samtools index -b -@32 $NodupX $baiNoX
#This bai files needs to be here or it will not run correctly in the haplotypecaller. 

#Here the GATK pipelinestarts and haplotypecaller is the only step needs to be in the loop. 
#Due to the rest being combined to one file in the step after this in hte genomicsDBimport step.
$gatk HaplotypeCaller -R $ref_bwa_mem -I $NodupX -O $Haplotype -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation

#Loop ends here due to the files being combined for the genmics db and will stay like that until at least the applyVQSR function.
#THe loop can be reintroduced later to 
rm $ALignedX
#rm $NodupX
done

# GenomicsDBImport approach
# Generate input_for_GenomicDBImport.txt, a tab-separated file that contains sample name and sample paths so the perl and perl2 script can work properly.
# input_for_GenomicDBImport.txt should look like this: This file is created in this format.

#GTEX-13PLJ /media/god/data1/Emil/haplotypecaller/GTEX-13PLJ-0003-SM-6WSCN.Haplotype.g.vcf.gz
#GTEX-UPIC /media/god/data1/Emil/haplotypecaller/GTEX-UPIC-0004-SM-5SOEF.Haplotype.g.vcf.gz
#GTEX-ZZPU /media/god/data1/Emil/haplotypecaller/GTEX-ZZPU-0003-SM-6WBUC.Haplotype.g.vcf.gz


 

# Also generate contigs.txt, a 1-column file with all the contigs from the g.vcf files. Open a terminal and run the command where you want the file generated.

# It can be generated using the script below:


# Now we create a GenomicDBImport path for joint genotyping. This Perl GenomicsBDImport script creates seperate databases per contig listed in the 'contigs.list' (here called 'contigs.txt'. We then genotype using GenotypeGCVFs across all of those databases to produce .vcf files for each contig. We then merge the vcfs into a final master file.

# Remember to make a folder called 'tmp'.

 

################################################# IMPORTANT ##########################################################################

# This is moved to perl.sh so the code will be found there with comments. For this part until Variantcaller Runs in two different scripts perlsh and perl2.sh . Also very important to fix the required txt files. 
#cut -f2 /media/god/data1/genome_indexes/bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict | cut -d':' -f2 | cut -d'_' -f1 | sort | uniq | grep -v 'Un' | grep -v 'EBV' | grep -v '1.6' > contigs.txt 
#run this in the terminal and make sure to not have a folder och file named chrX or any chr[int] names.
#contigs.txt should look like this:
#chr1
#chr10
#chr11
#chr12
#chr13
#chr14
# IF you are doing just 1 chromosome then you can do just a file with the chrX and that will be enough. 
#IMPORTANT is that the files you are working on only contain the X chromosome.
################################################# IMPORTANT ##########################################################################


# Once genotyped using this script we receive one vcf per contig and from there we need to collect all of those vcfs into a single file.

# Here I run a MergeVCF script, merging the output from the above GenotypeGVCFs script.

# first make a file with the location of all outputs from GenotypeGVCFs.

# input_variant_files.list looks like the following:
#Changed this to have a pointer to each file this is not optimised for a large number of samples then you need to use the list structure hinted at above.
#in3=/media/god/data1/Emil/haplotypecaller/GTEX-13PLJ-0003-SM-6WSCN.Haplotype.g.vcf.gz
#in2=/media/god/data1/Emil/haplotypecaller/GTEX-UPIC-0004-SM-5SOEF.Haplotype.g.vcf.gz
#in1=/media/god/data1/Emil/haplotypecaller/GTEX-ZZPU-0003-SM-6wBUC.Haplotype.g.vcf.gz
# Then run MergeVcfs to merge all the contig-split files into one vcf.
#$gatk MergeVcfs -I $in1 -I $in2 -I $in3 -O allSamples.vcf- This step is not done and due to not needing this due to only focusing on the X chromosome. 
# The merge part of the script is only needed when you have multiple chromosomes from the genomicsDb perl script.  


#Fix during friday to have the correct inputs so you can handle the data later. 
#$gatk GenotypeGVCFs -R $ref_bwa_mem -V Input.Xchr.g.vcf.gz -O $sample.GenotypeX.g.vcg.gz this part is done in the perl script check perl2 if you need more info on it.

#### VariantRecalibrator ####
# SNP modeling pass
out_vcf_GenotypeGVCFs=/media/god/data1/Emil/Control/fin/All_chrX.vcf.gz
variants_folder=/media/god/data1/genome_indexes/variants
vcf_temp=/media/god/data1/Emil/Control/tmp

output_vcf_SNP_recalfile_VariantRecalibrator=$vcf_temp/VarRecal.snps.recalfile.vcf
output_tranches_SNP_VariantRecalibrator=$vcf_temp/VarRecal.snps.VariantRecalibrator.tranches
output_rscript_SNP_VariantRecalibrator=$vcf_temp/VarRecal.snps.VariantRecalibrator.R

#$gatk VariantRecalibrator -R $ref_bwa_mem \
#-mode SNP \
#-AS \
#-an AS_QD -an AS_ReadPosRankSum -an AS_MQ -an AS_MQRankSum -an AS_SOR -an FS \
#--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $variants_folder/hapmap_3.3.hg38.vcf.gz \
#--resource:omni,known=false,training=true,truth=true,prior=12.0 $variants_folder/1000G_omni2.5.hg38.vcf.gz \
#--resource:1000G,known=false,training=true,truth=false,prior=10.0 $variants_folder/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
#--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $variants_folder/dbSNP_all_common_vcf/bkup/chr_00-common_all.vcf.gz \
#-tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.5 -tranche 99.3 -tranche 99.0 -tranche 90.0 \
#--variant $out_vcf_GenotypeGVCFs \
#-O $output_vcf_SNP_recalfile_VariantRecalibrator \
#--tranches-file $output_tranches_SNP_VariantRecalibrator \
#--rscript-file $output_rscript_SNP_VariantRecalibrator


# INDEL modeling pass
output_vcf_INDEL_recalfile_VariantRecalibrator=$vcf_temp/VarRecal.indel.recalfile.vcf
output_tranches_INDEL_VariantRecalibrator=$vcf_temp/VarRecal.indel.VariantRecalibrator.tranches
output_rscript_INDEL_VariantRecalibrator=$vcf_temp/VarRecal.indel.VariantRecalibrator.R

#$gatk VariantRecalibrator -R $ref_bwa_mem \
#-mode INDEL -AS \
#-an AS_QD -an AS_ReadPosRankSum -an AS_MQ -an AS_MQRankSum -an AS_SOR -an FS \
#--resource:1000G,known=false,training=true,truth=true,prior=12.0 $variants_folder/Mills_and_1000G_gold_standard.indels.hg38.vcf \
#--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $variants_folder/Homo_sapiens_assembly38.dbsnp138.vcf \
#-tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.5 -tranche 99.3 -tranche 99.0 -tranche 90.0 \
#--max-gaussians 4 \
#--variant $out_vcf_GenotypeGVCFs \
#-O $output_vcf_INDEL_recalfile_VariantRecalibrator \
#--tranches-file $output_tranches_INDEL_VariantRecalibrator \
#--rscript-file $output_rscript_INDEL_VariantRecalibrator



#### ApplyVSQR ####
# SNP filtering pass
output_vcf_SNP_ApplyVQSR=$vcf_temp/ApplyVSQRControl.snps.ApplyVQSR.vcf
#$gatk ApplyVQSR \
#-R $ref_bwa_mem \
#-mode SNP \
#--truth-sensitivity-filter-level 99.3 \
#-AS \
#--variant $out_vcf_GenotypeGVCFs \
#--recal-file $output_vcf_SNP_recalfile_VariantRecalibrator \
#--tranches-file $output_tranches_SNP_VariantRecalibrator \
#-O $output_vcf_SNP_ApplyVQSR


# INDEL filtering pass (use SNP filtering pass output as input)
output_vcf_INDEL_ApplyVQSR=$vcf_temp/ApplyVSQRControl.recalibrated.vcf

#$gatk ApplyVQSR -R $ref_bwa_mem \
#-mode INDEL \
#--truth-sensitivity-filter-level 99.3 \
#-AS \
#--variant $output_vcf_SNP_ApplyVQSR \
#--recal-file $output_vcf_INDEL_recalfile_VariantRecalibrator \
#--tranches-file $output_tranches_INDEL_VariantRecalibrator \
#-O $output_vcf_INDEL_ApplyVQSR


# To further analyse the output from the ApplyVQSR the use of the IGV software.

#$igv is the pointer for the folder of related igv data.

#This script is intended for launch on *nix machines

#-Xmx4g indicates 4 gb of memory, adjust number up or down as needed
#Add the flag -Ddevelopment = true to use features still in development
#Add the flag -Dsun.java2d.uiScale=2 for HiDPI displays
#prefix=`dirname $(readlink $0 || echo $0)`

# Check whether or not to use the bundled JDK
#written in command line is /media/god/data1/software/IGV_Linux_2.8.0/igv.sh
#To start the igv program.


# in here might be a place to put the Tiddid but tiddid is confined to complete genenome analysis so it cant be run on just the Xchr. 
# The tiddit is with the most probability optimised for the complete genome that will the cause issue when using the incomplete. 

# to do with tiddit first check if you can only run the samples of known non-mosaic and the after maybe run a handful of controls or if you should just leave that be and have it run all.
# also have to pull the tiddit from github to work and then install talk with Colm and Björn about this. 

#Run the tiddit for the 3 known and then on an extra 5 to have as control due to limited space on the harddrive.
