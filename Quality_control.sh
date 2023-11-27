#!/bin/bash


cd NanoporeFromBRC/paper3/Li_sequencing_sup

# QC for the R10.4.1 ligation sequencing with guppy -sup basecaller
export PYTHONPATH=/programs/nanofilt-2.8.0/lib/python3.9/site-packages
export PATH=/programs/nanofilt-2.8.0/bin:$PATH
# use three parallel screens to do the filtration is faster
#for i in *.fastq.gz; do gunzip -c $i | NanoFilt -q 10 -l 200 | gzip > ${i%%.*}_q10_len200.fastq.gz; done
gunzip -c Li_barcode01.fastq.gz | NanoFilt -q 10 -l 200 | gzip > Li_barcode01_q10_len200.fastq.gz
gunzip -c Li_barcode02.fastq.gz | NanoFilt -q 10 -l 200 | gzip > Li_barcode02_q10_len200.fastq.gz
gunzip -c Li_barcode03.fastq.gz | NanoFilt -q 10 -l 200 | gzip > Li_barcode03_q10_len200.fastq.gz

# visualized the quality of after QC
export PYTHONPATH=/programs/nanoplot-1.38.1/lib64/python3.6/site-packages:/programs/nanoplot-1.38.1/lib/python3.6/site-packages
export PATH=/programs/nanoplot-1.38.1/bin:$PATH
for i in *_q10_len200.fastq.gz; do echo NanoPlot -t 8 --fastq $i --plots dot --legacy hex -p ${i%_*}; done



# QC for the R9.4.1 ligation sequencing with guppy -sup basecaller
cd NanoporeFromBRC/paper3/BRC_sequencing_sup
# use parallel screens to do the filtration is faster
gunzip -c barcode02_passMinQ7.fastq.gz | NanoFilt -q 7 -l 200 | gzip > BRC_barcode02_q7_len200.fastq.gz
gunzip -c barcode03_passMinQ7.fastq.gz | NanoFilt -q 7 -l 200 | gzip > BRC_barcode03_q7_len200.fastq.gz
gunzip -c barcode04_passMinQ7.fastq.gz | NanoFilt -q 7 -l 200 | gzip > BRC_barcode04_q7_len200.fastq.gz
gunzip -c barcode05_passMinQ7.fastq.gz | NanoFilt -q 7 -l 200 | gzip > BRC_barcode05_q7_len200.fastq.gz
gunzip -c barcode07_passMinQ7.fastq.gz | NanoFilt -q 7 -l 200 | gzip > BRC_barcode07_q7_len200.fastq.gz
gunzip -c barcode08_passMinQ7.fastq.gz | NanoFilt -q 7 -l 200 | gzip > BRC_barcode08_q7_len200.fastq.gz
gunzip -c barcode09_passMinQ7.fastq.gz | NanoFilt -q 7 -l 200 | gzip > BRC_barcode09_q7_len200.fastq.gz

# visualized the quality of after QC
export PYTHONPATH=/programs/nanoplot-1.38.1/lib64/python3.6/site-packages:/programs/nanoplot-1.38.1/lib/python3.6/site-packages
export PATH=/programs/nanoplot-1.38.1/bin:$PATH
for i in *_q7_len200.fastq.gz; do echo NanoPlot -t 8 --fastq $i --plots dot --legacy hex -p ${i%_*}; done

# Assemble the sequences
cd /home/nw323/NanoporeFromBRC/paper3/
## Assemble the samples by lake - assemble R10.4.1 and R9.4.1 separatedly, and then merge
export PATH=/programs/Flye-2.9.2/bin:$PATH
### assemble with flye
flye --meta --nano-hq ./Li_sequencing_sup/Li_barcode02_q10_len200.fastq.gz ./Li_sequencing_sup/Li_barcode03_q10_len200.fastq.gz \
 --out-dir /home/nw323/NanoporeFromBRC/paper3/flye_Li_barcode0203 --threads 16
flye --meta --nano-hq ./BRC_sequencing_sup/BRC_barcode05_q7_len200.fastq.gz ./BRC_sequencing_sup/BRC_barcode09_q7_len200.fastq.gz \
--out-dir /home/nw323/NanoporeFromBRC/paper3/flye_BRC_barcode0507 --threads 16 
### first check draft bin when assembled separately
perl /programs/MaxBin-2.2.7/run_MaxBin.pl -contig flye_Li_barcode0203/assembly.fasta \
 -abund flye_Li_barcode0203/assembly_info.txt \
  -thread 16 -out MaxBin2_flye_Li_barcode0203 #out is output file head. 
perl /programs/MaxBin-2.2.7/run_MaxBin.pl -contig flye_BRC_barcode0507/assembly.fasta \
 -abund flye_BRC_barcode0507/assembly_info.txt \
  -thread 16 -out MaxBin2_flye_BRC_barcode0507 #out is output file head. 

mkdir R9_CN_coassemble
mv *BRC_barcode0507 R9_CN_coassemble

mkdir R10_CN_coassemble
mv *Li_barcode0203 R10_CN_coassemble

# use three rounds of racoon for polishing - Racon needs to be run on gen2
## use for loop to do racoon for each assembly
export PATH=/programs/minimap2-2.26:$PATH
export PATH=/programs/racon-1.5.0/bin:$PATH
allreads=R9_CN_AllReads.fastq.gz
assembly_file=R9_CN_coassemble/flye_BRC_barcode0507/assembly.fasta
cat ./BRC_sequencing_sup/BRC_barcode05_q7_len200.fastq.gz ./BRC_sequencing_sup/BRC_barcode09_q7_len200.fastq.gz \
> ${allreads}
for i in 1 2 3
do
    minimap2 -ax map-ont ${assembly_file} ${allreads} > ${allreads%_*}_aln_racon$[$i-1].sam
    racon -t 42 ${allreads} ${allreads%_*}_aln_racon$[$i-1].sam ${assembly_file} > racon${i}_${allreads%_*}.fasta
    assembly_file=racon${i}_${allreads%_*}.fasta
    echo "one round done"
done

allreads=R10_CN_AllReads.fastq.gz
assembly_file=R10_CN_coassemble/flye_Li_barcode0203/assembly.fasta
cat ./Li_sequencing_sup/Li_barcode02_q10_len200.fastq.gz ./Li_sequencing_sup/Li_barcode03_q10_len200.fastq.gz \
> ${allreads}
for i in 1 2 3
do
    minimap2 -ax map-ont ${assembly_file} ${allreads} > ${allreads%_*}_aln_racon$[$i-1].sam
    racon -t 42 ${allreads} ${allreads%_*}_aln_racon$[$i-1].sam ${assembly_file} > racon${i}_${allreads%_*}.fasta
    assembly_file=racon${i}_${allreads%_*}.fasta
    echo "one round done"
done

# use two rounds of medaka - need to use the GPU server and work in /workdir
##### for Guppy sup R9, the proper model is r941_min_sup_g507
##### for Guppy sup R10, the proper model is r1041_e82_400bps_sup_v4.2.0
mkdir /workdir/$USER/
mkdir /workdir/$USER/mydata
cp R9_CN_AllReads.fastq.gz /workdir/$USER/mydata/
cp racon3_R9_CN.fasta /workdir/$USER/mydata/
cp R10_CN_AllReads.fastq.gz /workdir/$USER/mydata/
cp racon3_R10_CN.fasta /workdir/$USER/mydata/
BASECALLS=R10_CN_AllReads.fastq.gz
DRAFT=racon3_R10_CN.fasta
OUTDIR=medaka1_R10_CN
model=r1041_e82_400bps_sup_v4.2.0
docker1 run -u root --rm -v /workdir/$USER/mydata/:/data -w /data  --gpus all ontresearch/medaka:latest medaka_consensus  -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR}  -m ${model}

BASECALLS=R10_CN_AllReads.fastq.gz
DRAFT=medaka1_R10_CN/consensus.fasta
OUTDIR=medaka2_R10_CN
model=r1041_e82_400bps_sup_v4.2.0
docker1 run -u root --rm -v /workdir/$USER/mydata/:/data -w /data  --gpus all ontresearch/medaka:latest medaka_consensus  -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR}  -m ${model}

# Maxbin2
READS=R10_CN_AllReads.fastq.gz
DRAFT=medaka2_R10_CN/consensus.fasta
PREFIX=maxbin_medaka_R10_CN
perl /programs/MaxBin-2.2.7/run_MaxBin.pl -contig ${DRAFT} \
 -reads ${READS} \
  -thread 16 -out ${PREFIX} #out is output file head. 

READS=R9_CN_AllReads.fastq.gz
DRAFT=medaka2_R9_CN/consensus.fasta
PREFIX=maxbin_medaka_R9_CN
perl /programs/MaxBin-2.2.7/run_MaxBin.pl -contig ${DRAFT} \
 -reads ${READS} \
  -thread 32 -out ${PREFIX} #out is output file head. 

# Maxbin2 on R9R10 contigs
cat S2_2_Medaka/medaka2_R10_CN/consensus.fasta S2_2_Medaka/medaka2_R9_CN/consensus.fasta > R9R10_CN_polished_contigs.fasta
DRAFT=R9R10_CN_polished_contigs.fasta
PREFIX=maxbin_medaka_R9R10_CN
perl /programs/MaxBin-2.2.7/run_MaxBin.pl -contig ${DRAFT} \
 -reads1 R9_CN_AllReads.fastq.gz -reads2 R10_CN_AllReads.fastq.gz \
  -thread 32 -out ${PREFIX} #out is output file head. 

----------------------- The files were moved to folders to better management ----------------
----------------------- Coaasemble output were put into S1_Flye, similar for Polish and Bin ----------------

# prepare the sorted bam file for 
READS=R9_CN_AllReads.fastq.gz
DRAFT=S2_2_Medaka/medaka2_R9_CN/consensus.fasta
bwa index ${DRAFT}
bwa mem -x ont2d -t 32 -M ${DRAFT} ${READS} > ${READS: 0 :5}_bwa.sam
samtools view -@ 32 -bS ${READS: 0 :5}_bwa.sam > ${READS: 0 :5}_bwa.bam
samtools sort -@ 32 ${READS: 0 :5}_bwa.bam > ${READS: 0 :5}_bwa.sort.bam
singularity run --bind $PWD --pwd $PWD /programs/metabat-2.16/metabat.sif runMetaBat.sh -s 500000 ${DRAFT} ${READS: 0 :5}_bwa.sort.bam 
mv consensus.fasta.depth.txt ${READS: 0 :5}_consensus.fasta.depth.txt
mv consensus.fasta.metabat-* ${READS: 0 :5}_metabat

READS=R10_CN_AllReads.fastq.gz
DRAFT=S2_2_Medaka/medaka2_R10_CN/consensus.fasta
bwa index ${DRAFT}
bwa mem -x ont2d -t 32 -M ${DRAFT} ${READS} > ${READS: 0 :6}_bwa.sam
samtools view -@ 32 -bS ${READS: 0 :6}_bwa.sam > ${READS: 0 :6}_bwa.bam
samtools sort -@ 32 ${READS: 0 :6}_bwa.bam > ${READS: 0 :6}_bwa.sort.bam
singularity run --bind $PWD --pwd $PWD /programs/metabat-2.16/metabat.sif runMetaBat.sh -s 500000 ${DRAFT} ${READS: 0 :6}_bwa.sort.bam
mv consensus.fasta.depth.txt ${READS: 0 :6}_consensus.fasta.depth.txt
mv consensus.fasta.metabat-* ${READS: 0 :6}_metabat

# # run vamb on R9R10 reads (it doesn't work for DAS)
# export PYTHONPATH=/programs/vamb-4.1.3/lib/python3.9/site-packages:/programs/vamb-4.1.3/lib64/python3.9/site-packages:/programs/vamb-4.1.3
# export PATH=/programs/vamb-4.1.3/bin:$PATH
# cp S2_2_Medaka/medaka2_R10_CN/consensus.fasta ./R10_CN_consensus.fasta
# cp S2_2_Medaka/medaka2_R9_CN/consensus.fasta ./R9_CN_consensus.fasta
# python concatenate.py R9R10_CN_polished_contigs_vamb.fasta R10_CN_consensus.fasta R9_CN_consensus.fasta --nozip
# cat R10_CN_AllReads.fastq.gz R9_CN_AllReads.fastq.gz > R9R10_CN_AllReads.fastq.gz
# DRAFT=R9R10_CN_polished_contigs_vamb.fasta
# bwa index ${DRAFT}
# bwa mem -x ont2d -t 32 -M ${DRAFT} R9R10_CN_AllReads.fastq.gz > ${DRAFT: 0 :8}_bwa.sam
# samtools view -@ 32 -bS ${DRAFT: 0 :8}_bwa.sam > ${DRAFT: 0 :8}_bwa.bam
# samtools sort -@ 32 ${DRAFT: 0 :8}_bwa.bam > ${DRAFT: 0 :8}_bwa.sort.bam

# CONTIG=R9R10_CN_polished_contigs_vamb.fasta
# BAMFILE=R9R10_CN_bwa.sort.bam  
# OUTPUT=R9R10_CN_Vamb
# vamb -o C --minfasta 500000 --outdir ${OUTPUT} --fasta ${CONTIG} --bamfiles ${BAMFILE}

# run vamb on R10 reads (R9 doesn't work because of small amount of contigs) 
CONTIG=R10_CN_consensus.fasta
BAMFILE=R10_CN_bwa.sort.bam  
OUTPUT=R10_CN_vamb
vamb --minfasta 500000 --outdir ${OUTPUT} --fasta ${CONTIG} --bamfiles ${BAMFILE}

# copy the files to the /workdir for DAStool
cd /home/nw323/NanoporeFromBRC/paper3/S3_1_Maxbin/maxbin_R9_CN
/home/nw323/NanoporeFromBRC/paper3/Fasta_to_Contig2Bin.sh -e fasta > maxbin_R9_CN_contigs2bin.tsv
sleep 1
cd /home/nw323/NanoporeFromBRC/paper3/S3_1_Maxbin/maxbin_R10_CN
/home/nw323/NanoporeFromBRC/paper3/Fasta_to_Contig2Bin.sh -e fasta > maxbin_R10_CN_contigs2bin.tsv
sleep 1
cd /home/nw323/NanoporeFromBRC/paper3/S3_2_metabat/R9_CN_metabat
/home/nw323/NanoporeFromBRC/paper3/Fasta_to_Contig2Bin.sh -e fa > metabat_R9_CN_contigs2bin.tsv
sleep 1
cd /home/nw323/NanoporeFromBRC/paper3/S3_2_metabat/R10_CN_metabat
/home/nw323/NanoporeFromBRC/paper3/Fasta_to_Contig2Bin.sh -e fa > metabat_R10_CN_contigs2bin.tsv
sleep 1
cd /home/nw323/NanoporeFromBRC/paper3/S3_3_Vamb/R10_CN_vamb/bins
/home/nw323/NanoporeFromBRC/paper3/Fasta_to_Contig2Bin.sh -e fna > vamb_R10_CN_contigs2bin.tsv
sleep 1

cd /workdir/nw323
singularity run -B /workdir/$USER --pwd /workdir/$USER /programs/DAS_Tool-1.1.6/das_tools.sif DAS_Tool \
-i metabat_R9_CN_contigs2bin.tsv,maxbin_R9_CN_contigs2bin.tsv \
-l metabat,maxbin \
-c R9_CN_consensus.fasta \
-o R9_CN_DAStool

singularity run -B /workdir/$USER --pwd /workdir/$USER /programs/DAS_Tool-1.1.6/das_tools.sif DAS_Tool \
-i metabat_R10_CN_contigs2bin.tsv,maxbin_R10_CN_contigs2bin.tsv,vamb_R10_CN_contigs2bin.tsv  \
-l metabat,maxbin,vamb \
-c R10_CN_consensus.fasta \
-o R10_CN_DAStool





