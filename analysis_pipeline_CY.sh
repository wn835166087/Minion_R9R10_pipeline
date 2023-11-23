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

# Maxbin2
READS=R10_CY_AllReads.fastq.gz
DRAFT=S2_2_Medaka/medaka2_R10_CY/consensus.fasta
PREFIX=maxbin_medaka_R10_CY
perl /programs/MaxBin-2.2.7/run_MaxBin.pl -contig ${DRAFT} \
 -reads ${READS} \
  -thread 32 -out ${PREFIX} #out is output file head. 

READS=R9_CY_AllReads.fastq.gz
DRAFT=S2_2_Medaka/medaka2_R9_CY/consensus.fasta
PREFIX=maxbin_medaka_R9_CY
perl /programs/MaxBin-2.2.7/run_MaxBin.pl -contig ${DRAFT} \
 -reads ${READS} \
  -thread 32 -out ${PREFIX} #out is output file head. 

# Metabat
# prepare the sorted bam file for 
READS=R9_CY_AllReads.fastq.gz
DRAFT=S2_2_Medaka/medaka2_R9_CY/consensus.fasta
bwa index ${DRAFT}
bwa mem -x ont2d -t 32 -M ${DRAFT} ${READS} > ${READS: 0 :5}_bwa.sam
samtools view -@ 32 -bS ${READS: 0 :5}_bwa.sam > ${READS: 0 :5}_bwa.bam
samtools sort -@ 32 ${READS: 0 :5}_bwa.bam > ${READS: 0 :5}_bwa.sort.bam
singularity run --bind $PWD --pwd $PWD /programs/metabat-2.16/metabat.sif runMetaBat.sh -s 500000 ${DRAFT} ${READS: 0 :5}_bwa.sort.bam 
mv consensus.fasta.depth.txt ${READS: 0 :5}_consensus.fasta.depth.txt
mv consensus.fasta.metabat-* ${READS: 0 :5}_metabat

READS=R10_CY_AllReads.fastq.gz
DRAFT=S2_2_Medaka/medaka2_R10_CY/consensus.fasta
bwa index ${DRAFT}
bwa mem -x ont2d -t 32 -M ${DRAFT} ${READS} > ${READS: 0 :6}_bwa.sam
samtools view -@ 32 -bS ${READS: 0 :6}_bwa.sam > ${READS: 0 :6}_bwa.bam
samtools sort -@ 32 ${READS: 0 :6}_bwa.bam > ${READS: 0 :6}_bwa.sort.bam
singularity run --bind $PWD --pwd $PWD /programs/metabat-2.16/metabat.sif runMetaBat.sh -s 500000 ${DRAFT} ${READS: 0 :6}_bwa.sort.bam
mv consensus.fasta.depth.txt ${READS: 0 :6}_consensus.fasta.depth.txt
mv consensus.fasta.metabat-* ${READS: 0 :6}_metabat

# Vamb
export PYTHONPATH=/programs/vamb-4.1.3/lib/python3.9/site-packages:/programs/vamb-4.1.3/lib64/python3.9/site-packages:/programs/vamb-4.1.3
export PATH=/programs/vamb-4.1.3/bin:$PATH
# R9 reads
CONTIG=R9_CY_consensus.fasta
BAMFILE=R9_CY_bwa.sort.bam  
OUTPUT=R9_CY_vamb
vamb --minfasta 500000 --outdir ${OUTPUT} --fasta ${CONTIG} --bamfiles ${BAMFILE}

# R10 reads
CONTIG=R10_CY_consensus.fasta
BAMFILE=R10_CY_bwa.sort.bam  
OUTPUT=R10_CY_vamb
vamb --minfasta 500000 --outdir ${OUTPUT} --fasta ${CONTIG} --bamfiles ${BAMFILE}

#
cd /home/nw323/NanoporeFromBRC/paper3/S3_1_Maxbin/maxbin_R9_CY
/home/nw323/NanoporeFromBRC/paper3/Fasta_to_Contig2Bin.sh -e fasta > maxbin_R9_CY_contigs2bin.tsv
sleep 1
cd /home/nw323/NanoporeFromBRC/paper3/S3_1_Maxbin/maxbin_R10_CY
/home/nw323/NanoporeFromBRC/paper3/Fasta_to_Contig2Bin.sh -e fasta > maxbin_R10_CY_contigs2bin.tsv
sleep 1
cd /home/nw323/NanoporeFromBRC/paper3/S3_2_metabat/R9_CY_metabat
/home/nw323/NanoporeFromBRC/paper3/Fasta_to_Contig2Bin.sh -e fa > metabat_R9_CY_contigs2bin.tsv
sleep 1
cd /home/nw323/NanoporeFromBRC/paper3/S3_2_metabat/R10_CY_metabat
/home/nw323/NanoporeFromBRC/paper3/Fasta_to_Contig2Bin.sh -e fa > metabat_R10_CY_contigs2bin.tsv
sleep 1
cd /home/nw323/NanoporeFromBRC/paper3/S3_3_Vamb/R10_CY_vamb/bins
/home/nw323/NanoporeFromBRC/paper3/Fasta_to_Contig2Bin.sh -e fna > vamb_R10_CY_contigs2bin.tsv
sleep 1
cd /home/nw323/NanoporeFromBRC/paper3/S3_3_Vamb/R9_CY_vamb/bins
/home/nw323/NanoporeFromBRC/paper3/Fasta_to_Contig2Bin.sh -e fna > vamb_R9_CY_contigs2bin.tsv
sleep 1

cd /workdir/nw323
singularity run -B /workdir/$USER --pwd /workdir/$USER /programs/DAS_Tool-1.1.6/das_tools.sif DAS_Tool \
-i metabat_R9_CY_contigs2bin.tsv,maxbin_R9_CY_contigs2bin.tsv \
-l metabat,maxbin \
-c R9_CY_consensus.fasta \
-o R9_CY_DAStool

singularity run -B /workdir/$USER --pwd /workdir/$USER /programs/DAS_Tool-1.1.6/das_tools.sif DAS_Tool \
-i metabat_R10_CY_contigs2bin.tsv,maxbin_R10_CY_contigs2bin.tsv,vamb_R10_CY_contigs2bin.tsv  \
-l metabat,maxbin,vamb \
-c R10_CY_consensus.fasta \
-o R10_CY_DAStool

