# merge the assembled contigs and rename it so that the later bining process is easier
export PYTHONPATH=/programs/vamb-4.1.3/lib/python3.9/site-packages:/programs/vamb-4.1.3/lib64/python3.9/site-packages:/programs/vamb-4.1.3
export PATH=/programs/vamb-4.1.3/bin:$PATH
python concatenate.py R9R10_CN_polished_contigs.fasta S2_2_Medaka/medaka2_R9_CY/consensus.fasta S2_2_Medaka/medaka2_R10_CY/consensus.fasta --nozip

# check if the R9 and R10 contigs has simialr contigs
export PATH=/programs/bbmap-39.03:$PATH
dedupe.sh in=R9R10_CY_polished_contigs.fasta out=R9R10_CY_polished_contigs_dedup_97.fasta outd=R9R10_CY_polished_contigs_duplicates_97.fa minidentity=97 t=32

# move CY reads together to facilitate downstream analysis
mkdir CY_reads
cp ../Li_sequencing_sup/Li_barcode01_q10_len200.fastq.gz CY_reads/
cp ../BRC_sequencing_sup/BRC_barcode02_q7_len200.fastq.gz CY_reads/
cp ../BRC_sequencing_sup/BRC_barcode03_q7_len200.fastq.gz CY_reads/
cp ../BRC_sequencing_sup/BRC_barcode04_q7_len200.fastq.gz CY_reads/
cp ../BRC_sequencing_sup/BRC_barcode07_q7_len200.fastq.gz CY_reads/
cp ../BRC_sequencing_sup/BRC_barcode08_q7_len200.fastq.gz CY_reads/

# maxbin
READS_ls=CY_reads_list_4maxbin.txt
DRAFT=R9R10_CN_polished_contigs.fasta
PREFIX=maxbin_CY
perl /programs/MaxBin-2.2.7/run_MaxBin.pl -contig ${DRAFT} \
-reads_list ${READS_ls} \
-thread 32 -out ${PREFIX} #out is output file head.

# prepare bam files for metabat and vamb
DRAFT=R9R10_CY_polished_contigs.fasta
bwa index ${DRAFT}
for fq_file in CY_reads/*.fastq.gz
do
    file_name=${fq_file##*/}
    bwa mem -x ont2d -t 32 -M ${DRAFT} ${fq_file} > ${file_name%_q*}_bwa.sam
    samtools view -@ 32 -bS ${file_name%_q*}_bwa.sam > ${file_name%_q*}_bwa.bam
    samtools sort -@ 32 ${file_name%_q*}_bwa.bam > ${file_name%_q*}_bwa.sort.bam
done

# Vamb
export PYTHONPATH=/programs/vamb-4.1.3/lib/python3.9/site-packages:/programs/vamb-4.1.3/lib64/python3.9/site-packages:/programs/vamb-4.1.3
export PATH=/programs/vamb-4.1.3/bin:$PATH
# R9 reads
CONTIG=R9R10_CY_polished_contigs.fasta
OUTPUT=S3_3_Vamb
vamb --minfasta 500000 --outdir ${OUTPUT} --fasta ${CONTIG} --bamfiles *_bwa.sort.bam 
cd ${OUTPUT}/bins
for i in `ls vae*`; do mv $i vamb_CY_${i%.*}.fasta; done
/home/nw323/NanoporeFromBRC/paper3/Fasta_to_Contig2Bin.sh -e fasta > vamb_CY_contigs2bin.tsv
cd ..


# Metabat, R9 and R10 run together 
DRAFT=R9R10_CY_polished_contigs.fasta
singularity run --bind $PWD --pwd $PWD /programs/metabat-2.16/metabat.sif runMetaBat.sh -s 500000 ${DRAFT} *_bwa.sort.bam  > metabat97_CY_log.txt 2>&1
mv ${DRAFT}.metabat-* S3_2_Metabat_R9R10
mv ${DRAFT}.depth.txt S3_2_Metabat_R9R10/
cd S3_2_Metabat_R9R10
for i in `ls bin*`; do mv $i metabat_CY_${i%.*}.fasta; done
/home/nw323/NanoporeFromBRC/paper3/Fasta_to_Contig2Bin.sh -e fasta > metabat_CY_R9R10_contigs2bin.tsv
cd ..


-------------
# post das_tool, use the in-house script: dastool_getbinfasta.py to get the bin fastas
python dastool_getbinfasta.py ${PWD}/

