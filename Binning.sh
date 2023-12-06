#!/bin/bash

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

# merge the R9 and R10 contigs
export PYTHONPATH=/programs/vamb-4.1.3/lib/python3.9/site-packages:/programs/vamb-4.1.3/lib64/python3.9/site-packages:/programs/vamb-4.1.3
export PATH=/programs/vamb-4.1.3/bin:$PATH
python concatenate.py R9R10_CY_polished_contigs.fasta S2_2_Medaka/medaka2_R9_CY/consensus.fasta S2_2_Medaka/medaka2_R10_CY/consensus.fasta --nozip

# check if the R9 contigs and R10 has simialr contigs
export PATH=/programs/bbmap-39.03:$PATH
dedupe.sh in=R9R10_CY_polished_contigs.fasta out=R9R10_CY_polished_contigs_dedup_97.fasta outd=R9R10_CY_polished_contigs_duplicates_97.fa minidentity=97 t=16
# then blast the R9R10_CY_polished_contigs_duplicates_97.fa
# there are 8 duplicate contigs, of which 6 are Microcystis contigs. 
# However, removing them doesn't improve the bin quality much, but certainly they should be removed in the finalized pipeline 


DRAFT=R9R10_CY_polished_contigs_dedup_97.fasta
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
CONTIG=R9R10_CY_polished_contigs_dedup_97.fasta
OUTPUT=S3_3_Vamb
vamb --minfasta 500000 --outdir ${OUTPUT} --fasta ${CONTIG} --bamfiles *_bwa.sort.bam 
cd ${OUTPUT}/bins
for i in `ls vae*`; do mv $i vamb_CY_${i%.*}.fasta; done
/home/nw323/NanoporeFromBRC/paper3/Fasta_to_Contig2Bin.sh -e fasta > vamb_CY_contigs2bin.tsv
cd ../..

DRAFT=R9R10_CY_polished_contigs_dedup_97.fasta
singularity run --bind $PWD --pwd $PWD /programs/metabat-2.16/metabat.sif runMetaBat.sh -s 500000 ${DRAFT} *_bwa.sort.bam  > metabat_log.txt 2>&1
mv ${DRAFT}.metabat-* S3_2_Metabat/
mv ${DRAFT}.depth.txt S3_2_Metabat/
cd S3_2_Metabat
for i in `ls bin*`; do mv $i metabat_CY_${i%.*}.fasta; done
/home/nw323/NanoporeFromBRC/paper3/Fasta_to_Contig2Bin.sh -e fasta > metabat_CY_contigs2bin.tsv
cp ${DRAFT}.depth.txt ${DRAFT}.depth.tsv
cd ..


awk -F'\t' '{
    for (i = 4; i <= NF; i += 2) {
        file = "Sample4metabet"(i-2)".txt";
        print $1, $i >> file;
    }
}' S3_2_Metabat/${DRAFT}.depth.tsv

for file in Sample4metabet*.txt; do
    # Extract sample name from the file
    sample_name=$(awk 'NR==1 {print $2}' "$file")

    # Rename the file with the sample name
    mv "$file" "${sample_name}.txt"

    # Remove the header from the file
    tail -n +2 "${sample_name}.txt" > "temp_${sample_name}.txt"
    mv "temp_${sample_name}.txt" "${sample_name}.txt"
done

# maxbin
ABUND_ls=CY_abund_list_4maxbin.txt
DRAFT=R9R10_CY_polished_contigs_dedup_97.fasta
PREFIX=maxbin_CY
perl /programs/MaxBin-2.2.7/run_MaxBin.pl -contig ${DRAFT} \
-abund_list ${ABUND_ls} \
-thread 32 -out ${PREFIX} #out is output file head.
mkdir S3_1_Maxbin
cp maxbin_CY.* S3_1_Maxbin/
cd S3_1_Maxbin
/home/nw323/NanoporeFromBRC/paper3/Fasta_to_Contig2Bin.sh -e fasta > maxbin_CY_contigs2bin.tsv
cd ..

rm -r /workdir/nw323
mkdir /workdir/nw323
cp S3_1_Maxbin/*_contigs2bin.tsv /workdir/nw323/
cp S3_2_Metabat/*_contigs2bin.tsv /workdir/nw323/
cp S3_3_Vamb/bins/*_contigs2bin.tsv /workdir/nw323/
cp R9R10_CY_polished_contigs_dedup_97.fasta /workdir/nw323/

cd /workdir/nw323
singularity run -B /workdir/$USER --pwd /workdir/$USER /programs/DAS_Tool-1.1.6/das_tools.sif DAS_Tool \
-i metabat_CY_contigs2bin.tsv,maxbin_CY_contigs2bin.tsv,vamb_CY_contigs2bin.tsv \
--score_threshold=0.1 \
-l metabat,maxbin,vamb \
-c R9R10_CY_polished_contigs_dedup_97.fasta \
-o S4_Dastool_CY

mkdir S4_Dastool
mv S4_Dastool_CY* S4_Dastool
cp -r S4_Dastool ~/NanoporeFromBRC/paper3/CY/

cd ~/NanoporeFromBRC/paper3/CY/S4_Dastool
/home/nw323/NanoporeFromBRC/paper3/dastool_getbinfasta.py ../R9R10_CY_polished_contigs_dedup_97.fasta S4_Dastool_CY_DASTool_contig2bin.tsv
