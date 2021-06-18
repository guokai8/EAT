for i in *_trimmed.fq.gz; 
do 
    bismark -L 32 -q --unmapped --ambiguous --phred33-quals --genome_folder /data/NGS/genomes/hisat2_annotation/Rat_EN/ -p 10 --multicore 5 $i ; 
done
for i in *.bam
do
    samtools sort -@ 40 -o ${i%%.bam}.sort.bam $i
done
for i in *.sort.bam
do
    samtools view -h $i >${i%%.bam}.sam
done
for i in *.sam; 
do 
    bismark_methylation_extractor -s --bedGraph --counts --buffer_size 50G --cytosine_report --genome_folder /data/NGS/genomes/hisat2_annotation/Rat_EN/ $i;rm -rf *sort.txt; 
done
