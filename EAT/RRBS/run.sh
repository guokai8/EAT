#!bin/bash
#QC filtering
for i in *_R1_001.fastq.gz; 
do 
    n=${i%%_R1_001.fastq.gz}; 
    java -jar /data/biotools/Trimmomatic/trimmomatic-0.36.jar PE -threads 50 -phred33 $i ${n}_R2_001.fastq.gz ${n}_R1_clean.fq.gz ${n}_R1_unpair.fq.gz ${n}_R2_clean.fq.gz ${n}_R2_unpair.fq.gz ILLUMINACLIP:/data/biotools/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 HEADCROP:13;
done
#alignment
for i in *_R1_clean.fq.gz 
do
 hisat2 --dta-cufflinks --rna-strandness RF -x /data/NGS/genomes/hisat2_annotation/Rat_EN/genome --known-splicesite-infile /data/NGS/genomes/hisat2_annotation/Rat_EN/genome.ss -p 40\
     -1 $i -2 ${i%%_R1_clean.fq.gz}_R2_clean.fq.gz|samtools sort -@ 20 -o ${i%%_R1_clean.fq.gz}.bam -
done
echo "All done"
###count calcuate
featureCounts -p -T 20 -s 2 -t exon -g gene_id -a /data/NGS/genomes/hisat2_annotation/Rat_EN/genes.gtf -o counts.txt *.bam

