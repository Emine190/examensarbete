cat /media/god/My Book Duo1/Emil/Control/sample | while read sample
do

echo $sample
control=/media/god/LaCie/GTEx_-50_females_whole_genome_seq/batch2
Input1=/media/god/My Book Duo1/Emil/Control/untrimmed/$sample'_'1.fastq.gz
Input2=/media/god/My Book Duo1/Emil/Control/untrimmed/$sample'_'2.fastq.gz

#Have samtools view function where the genome is and becomes a sorted bamfile.
samtools view -b -@32 $control/$sample.cram | samtools sort -o $sample.bam -@32 -n - 

samtools fastq -@32 -1 /media/god/My Book Duo1/Emil/Control/untrimmed/$sample'_'1.fastq.gz -2 /media/god/My Book Duo1/Emil/Control/untrimmed/$sample'_'2.fastq.gz  $sample.bam

rm $sample.bam

#removes the bam file intermediates so they dont take up unneccesary space on the disks.
# Here pointers for the reference genome and the pointers for the BWA and fastq.
ref_bwa_mem=/media/god/data1/genome_indexes/bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
trimmedin1=/media/god/My Book Duo1/Emil/Control/trimmed/trimmed.$sample'_'1.fastq.gz
trimmedin2=/media/god/My Book Duo1/Emil/Control/trimmed/trimmed.$sample'_'2.fastq.gz
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
$fastp -i /media/god/My Book Duo1/Emil/Control/untrimmed/$sample'_'1.fastq.gz \
-I /media/god/My Book Duo1/Emil/Control/untrimmed/$sample'_'2.fastq.gz \
-o /media/god/My Book Duo1/Emil/Control/trimmed/trimmed.$sample'_'1.fastq.gz \
-O /media/god/My Book Duo1/Emil/Control/trimmed/trimmed.$sample'_'2.fastq.gz

# trimming of the data to remove primers from the sequences.
rm $Input1
rm $Input2
#burrows-Wheeler Alignment tool to align the trimmed data from fastp. 
#the Functions will not work if you take out the X chromosome before the BWA it needs to be done after.

#Program: bwa (alignment via Burrows-Wheeler transformation)
#Version: 0.7.17-r1188

bwa mem -t32 -M -R "@RG\tID:$ID\tSM:$SM\tPL:ILLUMINA:150" $ref_bwa_mem\
/media/god/My Book Duo1/Emil/Control/trimmed/trimmed.$sample'_'1.fastq.gz \
/media/god/My Book Duo1/Emil/Control/trimmed/trimmed.$sample'_'1.fastq.gz \ | samtools sort -o /media/god/data1/Emil/Control/aligned/$sample.aligned.bam -@32 -

#Run an index to run the view functions to create just the X chromosome file.
#Creates the pointer for the index file for the index function so the bai files get to the bam files.
rm $trimmedin1
rm $trimmedin2

bai=/media/god'/My Book Duo1'/Emil/Control/aligned/$sample.aligned.bam.bai
#Rerun all the control samples from here in the control group. Something didnt go correctly from here. The pointer to bai was wrong.
samtools index -b -@32 /media/god/data1/Emil/Control/aligned/$sample.aligned.bam \
/media/god'/My Book Duo1'/Emil/Control/aligned/$sample.aligned.bam.bai
