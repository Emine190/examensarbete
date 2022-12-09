cat /media/MY/Files/Control/sample | while read sample
do
# This script is used to prep the samples for Tiddit so it just outputs bamfiles with their respective INdex file.
echo $sample
control=/media/files/Control/batch2
Input1=/media/MY/files/untrimmed/$sample'_'1.fastq.gz
Input2=/media/MY/files/untrimmed/$sample'_'2.fastq.gz

#Have samtools view function where the genome is and becomes a sorted bamfile.
samtools view -b -@32 $control/$sample.cram | samtools sort -o $sample.bam -@32 -n - 

samtools fastq -@32 -1 $Input1 -2 $Input2  $sample.bam

rm $sample.bam

#removes the bam file intermediates so they dont take up unneccesary space on the disks.
# Here pointers for the reference genome and the pointers for the BWA and fastq.
ref_bwa_mem=/media/genome_indexes/bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
trimmedin1=/media/MY/files/trimmed/trimmed.$sample'_'1.fastq.gz
trimmedin2=/media/MY/files/trimmed/trimmed.$sample'_'2.fastq.gz
ID=$(echo $sample | cut -d'-' -f1,2)
SM=$(echo $sample | cut -d'-' -f1,2)
Aligned=/media/MY/files/aligned/$sample.aligned.bam
echo $ID
echo $SM
echo $Input1
echo $Input2

# software
sambamba=/media/software/sambamba-0.7.1-linux-static
fastp=/media/software/fastp
gatk=/media/software/gatk-4.2.6.1/gatk
igv=/media/software/IGV_Linux_2.8.0

# quality trimming
$fastp -i $Input1 \
-I $Input2 \
-o $trimmedin1 \
-O $trimmedin2

# trimming of the data to remove primers from the sequences.
rm $Input1
rm $Input2
#burrows-Wheeler Alignment tool to align the trimmed data from fastp. 
#the Functions will not work if you take out the X chromosome before the BWA it needs to be done after.

#Program: bwa (alignment via Burrows-Wheeler transformation)
#Version: 0.7.17-r1188

bwa mem -t32 -M -R "@RG\tID:$ID\tSM:$SM\tPL:ILLUMINA:150" $ref_bwa_mem\
$trimmedin1 \
$trimmedin2 \ | samtools sort -o $Aligned -@32 -

#Run an index to run the view functions to create just the X chromosome file.
#Creates the pointer for the index file for the index function so the bai files get to the bam files.
rm $trimmedin1
rm $trimmedin2

bai=/media/MY/files/aligned/$sample.aligned.bam.bai
#Rerun all the control samples from here in the control group. Something didnt go correctly from here. The pointer to bai was wrong.
samtools index -b -@32 $Aligned \
/media/MY/files/aligned/$sample.aligned.bam.bai

done
