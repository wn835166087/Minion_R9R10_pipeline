# use the medium gen1: 24 cores, 128GB RAM
cd /home/nw323/NanoporeFromBRC/paper3/
cat ./BRC_sequencing_sup/BRC_barcode02_q7_len200.fastq.gz ./BRC_sequencing_sup/BRC_barcode03_q7_len200.fastq.gz \
./BRC_sequencing_sup/BRC_barcode04_q7_len200.fastq.gz ./BRC_sequencing_sup/BRC_barcode07_q7_len200.fastq.gz \
./BRC_sequencing_sup/BRC_barcode08_q7_len200.fastq.gz \
 > R9_CY_AllReads.fastq.gz

# Assemble the samples by lake - assemble R10.4.1 and R9.4.1 separatedly, and then merge
export PATH=/programs/Flye-2.9.2/bin:$PATH
flye --meta --nano-hq ./Li_sequencing_sup/Li_barcode01_q10_len200.fastq.gz \
 --out-dir ./flye_R10_CY --threads 16

flye --meta --nano-hq ./BRC_sequencing_sup/BRC_barcode02_q7_len200.fastq.gz ./BRC_sequencing_sup/BRC_barcode03_q7_len200.fastq.gz \
./BRC_sequencing_sup/BRC_barcode04_q7_len200.fastq.gz ./BRC_sequencing_sup/BRC_barcode07_q7_len200.fastq.gz \
./BRC_sequencing_sup/BRC_barcode08_q7_len200.fastq.gz \
--out-dir ./flye_R9_CY --threads 16 


# use for loop to do racoon for each assembly
# !!!racon need to be run on gen2 server!
export PATH=/programs/minimap2-2.26:$PATH
export PATH=/programs/racon-1.5.0/bin:$PATH

cp Li_sequencing_sup/Li_barcode01_q10_len200.fastq.gz R10_CY_AllReads.fastq.gz
allreads=R10_CY_AllReads.fastq.gz
assembly_file=flye_R10_CY/assembly.fasta
for i in 1 2 3
do
    minimap2 -ax map-ont ${assembly_file} ${allreads} > ${allreads%_*}_aln_racon$[$i-1].sam
    racon -t 48 ${allreads} ${allreads%_*}_aln_racon$[$i-1].sam ${assembly_file} > racon${i}_${allreads%_*}.fasta
    assembly_file=racon${i}_${allreads%_*}.fasta
    echo "one round done"
done

allreads=R9_CY_AllReads.fastq.gz
assembly_file=flye_R9_CY/assembly.fasta
for i in 1 2 3
do
    minimap2 -ax map-ont ${assembly_file} ${allreads} > ${allreads%_*}_aln_racon$[$i-1].sam
    racon -t 48 ${allreads} ${allreads%_*}_aln_racon$[$i-1].sam ${assembly_file} > racon${i}_${allreads%_*}.fasta
    assembly_file=racon${i}_${allreads%_*}.fasta
    echo "one round done"
done


# use two rounds of medaka - need to use the GPU server
##### for Guppy sup R9, the proper model is r941_min_sup_g507
##### for Guppy sup R10, the proper model is r1041_e82_400bps_sup_v4.2.0
mkdir /workdir/$USER/
mkdir /workdir/$USER/mydata
cp R9_CY_AllReads.fastq.gz /workdir/$USER/mydata/
cp racon3_R9_CY.fasta /workdir/$USER/mydata/
cp R10_CY_AllReads.fastq.gz /workdir/$USER/mydata/
cp racon3_R10_CY.fasta /workdir/$USER/mydata/

#docker1 run --rm -v /workdir/$USER/mydata/:/data -w /data  --gpus all ontresearch/medaka:latest medaka tools list\_models
BASECALLS=R9_CY_AllReads.fastq.gz
DRAFT=racon3_R9_CY.fasta
OUTDIR=medaka1_R9_CY
model=r941_min_sup_g507
docker1 run -u root --rm -v /workdir/$USER/mydata/:/data -w /data  --gpus all ontresearch/medaka:latest medaka_consensus  -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR}  -m ${model}

BASECALLS=R9_CY_AllReads.fastq.gz
DRAFT=medaka1_R9_CY/consensus.fasta
OUTDIR=medaka2_R9_CY
model=r941_min_sup_g507
docker1 run -u root --rm -v /workdir/$USER/mydata/:/data -w /data  --gpus all ontresearch/medaka:latest medaka_consensus  -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR}  -m ${model}



BASECALLS=R10_CY_AllReads.fastq.gz
DRAFT=racon3_R10_CY.fasta
OUTDIR=medaka1_R10_CY
model=r1041_e82_400bps_sup_v4.2.0
docker1 run -u root --rm -v /workdir/$USER/mydata/:/data -w /data  --gpus all ontresearch/medaka:latest medaka_consensus  -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR}  -m ${model}

BASECALLS=R10_CY_AllReads.fastq.gz
DRAFT=medaka1_R10_CY/consensus.fasta
OUTDIR=medaka2_R10_CY
model=r1041_e82_400bps_sup_v4.2.0
docker1 run -u root --rm -v /workdir/$USER/mydata/:/data -w /data  --gpus all ontresearch/medaka:latest medaka_consensus  -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR}  -m ${model}

---------------------- The files were moved to folders to better management ----------------
----------------------- Coaasemble output were put into S1_Flye, similar for Polish and Bin ----------------

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
