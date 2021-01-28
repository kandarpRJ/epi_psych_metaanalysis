library (exomePeak)

gtf<-c("hg38_gencode.v31.gtf")

######### human cerebrum

input <- c("PRJCA001180_human_cerebrum_1.input.bam", "PRJCA001180_human_cerebrum_2.input.bam")
ip <- c("PRJCA001180_human_cerebrum_1.m6A.bam", "PRJCA001180_human_cerebrum_2.m6A.bam")

result=exomepeak(GENE_ANNO_GTF=gtf, GENOME = "hg38", IP_BAM=ip, INPUT_BAM=input, SLIDING_STEP=50, OUTPUT_DIR="exomepeak_output", 
                 EXPERIMENT_NAME="human_cerebrum")

######### human cerebellum

input <- c("PRJCA001180_human_cerebellum_1.input.bam", "PRJCA001180_human_cerebellum_2.input.bam")
ip <- c("PRJCA001180_human_cerebellum_1.m6A.bam", "PRJCA001180_human_cerebellum_2.m6A.bam")

result=exomepeak(GENE_ANNO_GTF=gtf, GENOME = "hg38", IP_BAM=ip, INPUT_BAM=input, SLIDING_STEP=50, OUTPUT_DIR="exomepeak_output", 
                 EXPERIMENT_NAME="human_cerebellum")
