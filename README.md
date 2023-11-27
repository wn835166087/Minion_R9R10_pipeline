# MinION pipeline
This pipeline can be used to analyze the sequences generated via both R9 and R10 flowcells.
There are samples from two lakes, and the samples were co-assembled by two lakes. Here I only showed CN as an example
This pipeline can be used on Cornell Biohpc (https://biohpc.cornell.edu). The default server I used for this pipeline is medium gen1: 24 cores, 128GB RAM. However, there are a few softwares, specificly, racon and medaka need to be run specific servers (see instructions in analysis_pipeline_CN.txt). Therefore, instead of providing a shell script, here I only provide a text file to record the pipeline I used for analyzing the MinION data.
