cd NanoporeFromBRC/paper3/Li_sequencing_sup

# QC for the R10.4.1 ligation sequencing
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



# QC for the R9.4.1 ligation sequencing
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



# Assemble the samples by lake - only R10.4.1 - not necessary
# export PATH=/programs/Flye-2.9.2/bin:$PATH
# flye --meta --nano-hq Li_barcode01_q10_len200.fastq.gz --out-dir /home/nw323/NanoporeFromBRC/paper3/flye_Li_barcode01 --threads 16
# flye --meta --nano-hq Li_barcode02_q10_len200.fastq.gz Li_barcode03_q10_len200.fastq.gz --out-dir /home/nw323/NanoporeFromBRC/paper3/flye_Li_barcode0203 --threads 16

cd /home/nw323/NanoporeFromBRC/paper3/
# Assemble the samples by lake - R9.4.1 + R.10.4.1
export PATH=/programs/Flye-2.9.2/bin:$PATH
# use CN assembly to test if there's improvement when co-assemble with R10.4.1
flye --meta --nano-hq ./Li_sequencing_sup/Li_barcode02_q10_len200.fastq.gz ./Li_sequencing_sup/Li_barcode03_q10_len200.fastq.gz \
./BRC_sequencing_sup/BRC_barcode05_q7_len200.fastq.gz ./BRC_sequencing_sup/BRC_barcode09_q7_len200.fastq.gz \
--out-dir /home/nw323/NanoporeFromBRC/paper3/flye_LiBRC_CN --threads 16 
# draft bin it, check bin quality and compare with previous assembly
perl /programs/MaxBin-2.2.7/run_MaxBin.pl -contig flye_LiBRC_CN/assembly.fasta \
 -abund flye_LiBRC_CN/assembly_info.txt \
  -thread 16 -out MaxBin2_flye_LiBRC_CN #out is output file head. 


flye --meta --nano-hq ./Li_sequencing_sup/Li_barcode01_q10_len200.fastq.gz ./BRC_sequencing_sup/BRC_barcode02_q7_len200.fastq.gz \
./BRC_sequencing_sup/BRC_barcode03_q7_len200.fastq.gz ./BRC_sequencing_sup/BRC_barcode04_q7_len200.fastq.gz \
./BRC_sequencing_sup/BRC_barcode07_q7_len200.fastq.gz ./BRC_sequencing_sup/BRC_barcode08_q7_len200.fastq.gz \
--out-dir /home/nw323/NanoporeFromBRC/paper3/flye_LiBRC_CY --threads 16





# use two rounds of racoon for polishing
export PATH=/programs/racon-1.5.0/bin:$PATH
