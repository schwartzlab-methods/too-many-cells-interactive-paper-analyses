#!/bin/bash

DATA_DIR=./data/raw

### Prostate cancer
# LNCaP DMSO (control): GSM5155457
mkdir -p "$DATA_DIR/lncap_control"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5155nnn/GSM5155455/suppl/GSM5155455_LNCaP-DMSO_barcodes.tsv.gz -O "$DATA_DIR/lncap_control/barcodes.tsv.gz"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5155nnn/GSM5155455/suppl/GSM5155455_LNCaP-DMSO_features.tsv.gz -O "$DATA_DIR/lncap_control/features.tsv.gz"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5155nnn/GSM5155455/suppl/GSM5155455_LNCaP-DMSO_matrix.mtx.gz -O "$DATA_DIR/lncap_control/matrix.mtx.gz"
# LNCaP enzalutamide (anti-androgen): GSM5155457
mkdir "$DATA_DIR/lncap"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5155nnn/GSM5155457/suppl/GSM5155457_LNCaP-RESA_barcodes.tsv.gz -O "$DATA_DIR/lncap/barcodes.tsv.gz"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5155nnn/GSM5155457/suppl/GSM5155457_LNCaP-RESA_features.tsv.gz -O "$DATA_DIR/lncap/features.tsv.gz"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5155nnn/GSM5155457/suppl/GSM5155457_LNCaP-RESA_matrix.mtx.gz -O "$DATA_DIR/lncap/matrix.mtx.gz"

### T-ALL
# DND-41 DMSO (control): GSM4121361
mkdir "$DATA_DIR/dnd41_control"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4121nnn/GSM4121361/suppl/GSM4121361_DND41_DMSO_scRNA_barcodes.tsv.gz -O "$DATA_DIR/dnd41_control/barcodes.tsv.gz"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4121nnn/GSM4121361/suppl/GSM4121361_DND41_DMSO_scRNA_genes.tsv.gz -O "$DATA_DIR/dnd41_control/genes.tsv.gz"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4121nnn/GSM4121361/suppl/GSM4121361_DND41_DMSO_scRNA_matrix.mtx.gz -O "$DATA_DIR/dnd41_control/matrix.mtx.gz"
gunzip $DATA_DIR/dnd41_control/*
# DND-41 gamma secretase inhibitor (GSI) sustained: GSM4121364
mkdir "$DATA_DIR/dnd41_long"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4121nnn/GSM4121364/suppl/GSM4121364_DND41_GSI_Res_sustained_barcodes.tsv.gz -O "$DATA_DIR/dnd41_long/barcodes.tsv.gz"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4121nnn/GSM4121364/suppl/GSM4121364_DND41_GSI_Res_sustained_genes.tsv.gz -O "$DATA_DIR/dnd41_long/genes.tsv.gz"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4121nnn/GSM4121364/suppl/GSM4121364_DND41_GSI_Res_sustained_matrix.mtx.gz -O "$DATA_DIR/dnd41_long/matrix.mtx.gz"
gunzip $DATA_DIR/dnd41_long/*
# DND-41 gamma secretase inhibitor (GSI) sustained: GSM4121364
mkdir "$DATA_DIR/dnd41_short"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4121nnn/GSM4121362/suppl/GSM4121362_DND41_GSI_scRNA_barcodes.tsv.gz -O "$DATA_DIR/dnd41_short/barcodes.tsv.gz"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4121nnn/GSM4121362/suppl/GSM4121362_DND41_GSI_scRNA_genes.tsv.gz -O "$DATA_DIR/dnd41_short/genes.tsv.gz"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4121nnn/GSM4121362/suppl/GSM4121362_DND41_GSI_scRNA_matrix.mtx.gz -O "$DATA_DIR/dnd41_short/matrix.mtx.gz"
gunzip $DATA_DIR/dnd41_short/*

### Melanoma
# SK-MEL-28 DMSO (control): GSM4932163
mkdir "$DATA_DIR/skmel28_control"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4932nnn/GSM4932163/suppl/GSM4932163_sample5_barcodes.tsv.gz -O "$DATA_DIR/skmel28_control/barcodes.tsv.gz"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4932nnn/GSM4932163/suppl/GSM4932163_sample5_genes.tsv.gz -O "$DATA_DIR/skmel28_control/genes.tsv.gz"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4932nnn/GSM4932163/suppl/GSM4932163_sample5_matrix.mtx.gz -O "$DATA_DIR/skmel28_control/matrix.mtx.gz"
gunzip $DATA_DIR/skmel28_control/*
# SK-MEL-28 dabrafenin-treated 72h: GSM4932166
mkdir "$DATA_DIR/skmel28"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4932nnn/GSM4932166/suppl/GSM4932166_sample8_barcodes.tsv.gz -O "$DATA_DIR/skmel28/barcodes.tsv.gz"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4932nnn/GSM4932166/suppl/GSM4932166_sample8_genes.tsv.gz -O "$DATA_DIR/skmel28/genes.tsv.gz"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4932nnn/GSM4932166/suppl/GSM4932166_sample8_matrix.mtx.gz -O "$DATA_DIR/skmel28/matrix.mtx.gz"
gunzip $DATA_DIR/skmel28/*

### Breast cancer
# MDA-MB-231 epithelial adenocarcinoma DMSO (control): GSM4684556
mkdir "$DATA_DIR/mdamb231_control"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4684nnn/GSM4684556/suppl/GSM4684556_t0_barcodes.tsv.gz -O "$DATA_DIR/mdamb231_control/barcodes.tsv.gz"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4684nnn/GSM4684556/suppl/GSM4684556_t0_features.tsv.gz -O "$DATA_DIR/mdamb231_control/features.tsv.gz"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4684nnn/GSM4684556/suppl/GSM4684556_t0_matrix.mtx.gz -O "$DATA_DIR/mdamb231_control/matrix.mtx.gz"
# MDA-MB-231 doxorubicin treated 7-weeks: GSM4684557
mkdir "$DATA_DIR/mdamb231"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4684nnn/GSM4684557/suppl/GSM4684557_t7_barcodes.tsv.gz -O "$DATA_DIR/mdamb231/barcodes.tsv.gz"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4684nnn/GSM4684557/suppl/GSM4684557_t7_features.tsv.gz -O "$DATA_DIR/mdamb231/features.tsv.gz"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4684nnn/GSM4684557/suppl/GSM4684557_t7_matrix.mtx.gz -O "$DATA_DIR/mdamb231/matrix.mtx.gz"

### Lung cancer
# PC9 untreated: GSM3972651
mkdir "$DATA_DIR/pc9_control"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3972nnn/GSM3972651/suppl/GSM3972651_PC9D0_filtered_feature_bc_matrices.tar.gz -O "$DATA_DIR/pc9_control.tar.gz"
tar -xvf $DATA_DIR/pc9_control.tar.gz -C $DATA_DIR/pc9_control
mv $DATA_DIR/pc9_control/**/*.gz $DATA_DIR/pc9_control
rm -r $DATA_DIR/pc9_control/home
rm $DATA_DIR/pc9_control.tar.gz
# PC9 erlotinib (EGFR inhibitor) treated 3 days: GSM3972652
mkdir "$DATA_DIR/pc9"
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3972nnn/GSM3972652/suppl/GSM3972652_PC9D3Erl_filtered_feature_bc_matrices.tar.gz -O "$DATA_DIR/pc9.tar.gz"
tar -xvf $DATA_DIR/pc9.tar.gz -C $DATA_DIR/pc9
mv $DATA_DIR/pc9/**/*.gz $DATA_DIR/pc9
rm -r $DATA_DIR/pc9/home
rm $DATA_DIR/pc9.tar.gz
