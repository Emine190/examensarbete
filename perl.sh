#!/usr/bin/perl

################################################# IMPORTANT ##########################################################################

# These two script needs to be run separately from the other scripts (as it uses perl). So copy the below two scripts into a new .sh file, making sure that the file begins with #!/usr/bin/perl .step 1.

################################################# IMPORTANT ##########################################################################





use warnings;

use strict;

my $gatkcommand = "/media/god/data1/software/gatk-4.2.6.1/gatk --java-options ' -Djava.io.tmpdir=/tmp' GenomicsDBImport -R /media/god/data1/genome_indexes/bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --sample-name-map input_for_GenomicDBImport.txt --reader-threads 32 ";

my $contigfile = 'contigs.txt';

open my $contig, $contigfile or die "Could not open $contigfile: $!";

while( my $line = <$contig>)  {

chomp($line);

my $gatkfinalcommand = $gatkcommand."--genomicsdb-workspace-path ".$line." -L ".$line;

system($gatkfinalcommand);

}

 
# The raw .g.vcf file is then run through GenotypeGVCFs which performs joint genotyping on all samples together. step 2
####################################################################################################################
#!/usr/bin/perl

#use warnings;

#use strict;

#my $gatkcommand = "/media/god/data1/software/gatk-4.2.6.1/gatk --java-options '-Xmx64g -Djava.io.tmpdir=/tmp' GenotypeGVCFs -R /media/god/data1/genome_indexes/bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -G StandardAnnotation -G AS_StandardAnnotation";

#my $contigfile = 'contigs.txt';

#open my $contig, $contigfile or die "Could not open $contigfile: $!";

#while( my $line = <$contig>)  {

 #   chomp($line);

  #      my $gatkfinalcommand = $gatkcommand." --variant gendb://".$line." -O ./fin/All_".$line.".vcf.gz";

   #     system($gatkfinalcommand);

#}


 
