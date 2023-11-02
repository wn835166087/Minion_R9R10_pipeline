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





