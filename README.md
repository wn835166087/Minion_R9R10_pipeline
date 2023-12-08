# Analysis pipeline for MinION data sequenced by R9 and R10 flowcells.
This pipeline can be used to analyze the sequences generated via both R9 and R10 flowcells.

There are samples from two lakes, and the samples were co-assembled by two lakes. The pipeline are exactly the same for both lakes.

This pipeline can be used on Cornell Biohpc (https://biohpc.cornell.edu). The default server I used for this pipeline is medium gen1: 24 cores, 128GB RAM. The order of using the pipeline is Quality_Control -> Assemble_Polish -> Binning. For Quality_Control and Binning, the .sh file can be submitted to the server. However, the softwares for the contigs polish, specificly, racon and medaka need to be run specific servers. Therefore, instead of providing a shell script, here I only provide a text file to record the assemble and polish pipeline I used for analyzing the MinION data.

### Some other thought about the pipeline:
Note that for metabat and maxbin, one strategy is to lower the percentage identity when calculating the contig depth, then combine the depth information for binning. It helps with getting more high-quality MAGs. But doesn’t help much Microcystis MAGs. Therefore, It was not applied.

