#!/bin/bash
export JAVA_HOME=/home/tooba/.sdkman/candidates/java/current
export JAVA_CMD=$JAVA_HOME/bin/java
export NXF_VER=25.04.0
export PATH=$JAVA_HOME/bin:/home/tooba/bin:$PATH

source ~/miniconda3/etc/profile.d/conda.sh
conda activate r_env

/home/tooba/bin/nextflow run /home/tooba/Skin_Cell_Model/l-ctpred/PredictDb-nextflow/main.nf \
--gene_annotation /home/tooba/Skin_Cell_Model/l-ctpred/gencode.v49.annotation.gtf \
--snp_annotation /home/tooba/Skin_Cell_Model/l-ctpred/all_chrs.snp_annot_v1.txt \
--genotype /home/tooba/Skin_Cell_Model/l-ctpred/all_chrs.geno_v1.txt \
--gene_exp /home/tooba/Skin_Cell_Model/predictions/Venous_EC_predictions.csv \
--outdir /home/tooba/Skin_Cell_Model/predictdb_output/Venous_EC \
--nfolds 3 \
--max_cpus 30


