

export PYTHONPATH=/programs/vamb-4.1.3/lib/python3.9/site-packages:/programs/vamb-4.1.3/lib64/python3.9/site-packages:/programs/vamb-4.1.3
export PATH=/programs/vamb-4.1.3/bin:$PATH
python concatenate.py R9R10_CN_polished_contigs.fasta S2_2_Medaka/medaka2_R9_CY/consensus.fasta S2_2_Medaka/medaka2_R10_CY/consensus.fasta --nozip

# maxbin
READS_ls=CY_reads_list_4maxbin.txt
DRAFT=R9R10_CN_polished_contigs.fasta
PREFIX=maxbin_CY
perl /programs/MaxBin-2.2.7/run_MaxBin.pl -contig ${DRAFT} \
-reads_list ${READS_ls} \
-thread 32 -out ${PREFIX} #out is output file head.

# prepare bam files for metabat and vamb
DRAFT=R9R10_CN_polished_contigs.fasta
bwa index ${DRAFT}
for fq_file in BRC_sequencing_sup/*_q7_len200.fastq.gz
do
    file_name=${fq_file#*/}
    bwa mem -x ont2d -t 32 -M ${DRAFT} ${fq_file} > ${file_name%_q*}_bwa.sam
    samtools view -@ 32 -bS ${file_name%_q*}_bwa.sam > ${file_name%_q*}_bwa.bam
    samtools sort -@ 32 ${file_name%_q*}_bwa.bam > ${file_name%_q*}_bwa.sort.bam
done

DRAFT=R9R10_CN_polished_contigs.fasta
for fq_file in Li_sequencing_sup/*_q10_len200.fastq.gz
do
    file_name=${fq_file#*/}
    bwa mem -x ont2d -t 32 -M ${DRAFT} ${fq_file} > ${file_name%_q*}_bwa.sam
    samtools view -@ 32 -bS ${file_name%_q*}_bwa.sam > ${file_name%_q*}_bwa.bam
    samtools sort -@ 32 ${file_name%_q*}_bwa.bam > ${file_name%_q*}_bwa.sort.bam
done



-------------
# post das_tool, use the in-house script: dastool_getbinfasta.py to get the bin fastas
python dastool_getbinfasta.py ${PWD}/

