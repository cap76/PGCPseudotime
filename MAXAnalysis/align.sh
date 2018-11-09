bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Sequence/Bowtie2Index/genome -U /mnt/scratch/surani/cap76/MAX/SRR643235_trimmed.fq.gz -S /mnt/scratch/surani/cap76/MAX/SRR643235_trimmed.sam

samtools view -bS  /mnt/scratch/surani/cap76/MAX/SRR643235_trimmed.sam | samtools sort - /mnt/scratch/surani/cap76/MAX/SRR643235_trimmed_sorted
samtools index /mnt/scratch/surani/cap76/MAX/SRR643235_trimmed_sorted.bam /mnt/scratch/surani/cap76/MAX/SRR643235_trimmed_sorted.bai

bamCoverage -b  /mnt/scratch/surani/cap76/MAX/SRR643235_trimmed_sorted.bam  --binSize 20 --extendReads 150 --normalizeUsingRPKM -o  /mnt/scratch/surani/cap76/MAX/SRR643235_trimmed_sorted.bw
macs2 callpeak -t /mnt/scratch/surani/cap76/MAX/SRR643235_trimmed.sam -g hs -n  /mnt/scratch/surani/cap76/MAX/SRR643235_trimmed_sorted -q 0.05

