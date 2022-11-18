
cat /media/MY/files/Control/sample | while read sample
do

echo $sample
ref_bwa_mem=/media/files/genome_indexes/bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
gatk=/media/software/gatk-4.2.6.1/gatk

Haplotype=/media/MY/files/haplotypecaller/$sample.Haplotype.g.vcf.gz
#rm $Aligned
#before this step check the head and tail of the files so they are OK other wise redo.
baiNoX=/media/god'/My Book Duo1'/Emil/Control/NodupXchr/$sample.NodupXchr.bam.bai
samtools index -b -@32 /media/MY/files/NodupXchr/$sample.NodupXchr.bam \
/media/MY/files/NodupXchr/$sample.NodupXchr.bam.bai
#This bai files needs to be here or it will not run correctly in the haplotypecaller. 

#Here the GATK pipelinestarts and haplotypecaller is the only step needs to be in the loop. 
#Due to the rest being combined to one file in the step after this in hte genomicsDBimport step.
$gatk HaplotypeCaller -R $ref_bwa_mem -I /media/MY/files/NodupXchr/$sample.NodupXchr.bam \
-O /media/MY/files/Control/haplotypecaller/$sample.Haplotype.g.vcf.gz \
-ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation

#Loop ends here due to the files being combined for the genmics db and will stay like that until at least the applyVQSR function.
#THe loop can be reintroduced later to 
rm $ALignedX
#rm $NodupX
done

# GenomicsDBImport approach
# Generate input_for_GenomicDBImport.txt, a tab-separated file that contains sample name and sample paths so the perl and perl2 script can work properly.
# input_for_GenomicDBImport.txt should look like this: This file is created in this format.

#Woman1 /media/MY/files/haplotypecaller/Woman1.g.vcf.gz
#Woman2 /media/MY/files/haplotypecaller/Woman2.g.vcf.gz
#Woman3 /media/MY/files/haplotypecaller/Woman3.g.vcf.gz


 

# Also generate contigs.txt, a 1-column file with all the contigs from the g.vcf files. Open a terminal and run the command where you want the file generated.

# It can be generated using the script below:


# Now we create a GenomicDBImport path for joint genotyping. This Perl GenomicsBDImport script creates seperate databases per contig listed in the 'contigs.list' (here called 'contigs.txt'. We then genotype using GenotypeGCVFs across all of those databases to produce .vcf files for each contig. We then merge the vcfs into a final master file.

# Remember to make a folder called 'tmp'.

 

################################################# IMPORTANT ##########################################################################

# This is moved to perl.sh so the code will be found there with comments. For this part until Variantcaller Runs in two different scripts perlsh and perl2.sh . Also very important to fix the required txt files. 
#cut -f2 /media/files/genome_indexes/bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict | cut -d':' -f2 | cut -d'_' -f1 | sort | uniq | grep -v 'Un' | grep -v 'EBV' | grep -v '1.6' > contigs.txt 
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
#### VariantRecalibrator ####
# SNP modeling pass
out_vcf_GenotypeGVCFs=/media/MY/files/Control/fin/All_chrX.vcf.gz
variants_folder=/media/files/genome_indexes/variants
vcf_temp=/media/MY/files/Control/tmp

output_vcf_SNP_recalfile_VariantRecalibrator=$vcf_temp/VarRecal.snps.recalfile.vcf
output_tranches_SNP_VariantRecalibrator=$vcf_temp/VarRecal.snps.VariantRecalibrator.tranches
output_rscript_SNP_VariantRecalibrator=$vcf_temp/VarRecal.snps.VariantRecalibrator.R

$gatk VariantRecalibrator -R $ref_bwa_mem \
-mode SNP \
-AS \
-an AS_QD -an AS_ReadPosRankSum -an AS_MQ -an AS_MQRankSum -an AS_SOR -an FS \
--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $variants_folder/hapmap_3.3.hg38.vcf.gz \
--resource:omni,known=false,training=true,truth=true,prior=12.0 $variants_folder/1000G_omni2.5.hg38.vcf.gz \
--resource:1000G,known=false,training=true,truth=false,prior=10.0 $variants_folder/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $variants_folder/dbSNP_all_common_vcf/bkup/chr_00-common_all.vcf.gz \
-tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.5 -tranche 99.3 -tranche 99.0 -tranche 90.0 \
--variant $out_vcf_GenotypeGVCFs \
-O $output_vcf_SNP_recalfile_VariantRecalibrator \
--tranches-file $output_tranches_SNP_VariantRecalibrator \
--rscript-file $output_rscript_SNP_VariantRecalibrator


# INDEL modeling pass
output_vcf_INDEL_recalfile_VariantRecalibrator=$vcf_temp/VarRecal.indel.recalfile.vcf
output_tranches_INDEL_VariantRecalibrator=$vcf_temp/VarRecal.indel.VariantRecalibrator.tranches
output_rscript_INDEL_VariantRecalibrator=$vcf_temp/VarRecal.indel.VariantRecalibrator.R

$gatk VariantRecalibrator -R $ref_bwa_mem \
-mode INDEL -AS \
-an AS_QD -an AS_ReadPosRankSum -an AS_MQ -an AS_MQRankSum -an AS_SOR -an FS \
--resource:1000G,known=false,training=true,truth=true,prior=12.0 $variants_folder/Mills_and_1000G_gold_standard.indels.hg38.vcf \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $variants_folder/Homo_sapiens_assembly38.dbsnp138.vcf \
-tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.5 -tranche 99.3 -tranche 99.0 -tranche 90.0 \
--max-gaussians 4 \
--variant $out_vcf_GenotypeGVCFs \
-O $output_vcf_INDEL_recalfile_VariantRecalibrator \
--tranches-file $output_tranches_INDEL_VariantRecalibrator \
--rscript-file $output_rscript_INDEL_VariantRecalibrator



#### ApplyVSQR ####
# SNP filtering pass
output_vcf_SNP_ApplyVQSR=$vcf_temp/ApplyVSQRControl.snps.ApplyVQSR.vcf
$gatk ApplyVQSR \
-R $ref_bwa_mem \
-mode SNP \
--truth-sensitivity-filter-level 99.3 \
-AS \
--variant $out_vcf_GenotypeGVCFs \
--recal-file $output_vcf_SNP_recalfile_VariantRecalibrator \
--tranches-file $output_tranches_SNP_VariantRecalibrator \
-O $output_vcf_SNP_ApplyVQSR


# INDEL filtering pass (use SNP filtering pass output as input)
output_vcf_INDEL_ApplyVQSR=$vcf_temp/ApplyVSQRControl.recalibrated.vcf

$gatk ApplyVQSR -R $ref_bwa_mem \
-mode INDEL \
--truth-sensitivity-filter-level 99.3 \
-AS \
--variant $output_vcf_SNP_ApplyVQSR \
--recal-file $output_vcf_INDEL_recalfile_VariantRecalibrator \
--tranches-file $output_tranches_INDEL_VariantRecalibrator \
-O $output_vcf_INDEL_ApplyVQSR
