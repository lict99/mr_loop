#! /bin/bash

OUTPUT_DIR="./results/03"
RAW_BFILE="./data/SNParray/all/raw_data/22B0218D_23B0203E_23P1212B00A_24B0516F"
IMP_BFILE="./data/SNParray/all/impute_qc_data/final"
SEX_FILE="./results/02/sex_info.txt"

# Setting directory for the output files
if [ ! -d $OUTPUT_DIR ]; then
    mkdir -p $OUTPUT_DIR
fi

# Updating and filtering
plink \
    --bfile $RAW_BFILE \
    --update-sex $SEX_FILE \
    --mind 0.05 \
    --geno 0.05 \
    --make-bed \
    --out "$OUTPUT_DIR/data_sex_check"

# Sex check
plink \
    --bfile "$OUTPUT_DIR/data_sex_check" \
    --check-sex \
    --out "$OUTPUT_DIR/sex_check"

# Extraction of problematic samples
grep "PROBLEM" "$OUTPUT_DIR/sex_check.sexcheck" | awk '{print$1,$2}' >"$OUTPUT_DIR/sex_discrepancy.txt"

# Removing problematic samples
plink \
    --bfile $IMP_BFILE \
    --remove "$OUTPUT_DIR/sex_discrepancy.txt" \
    --make-bed \
    --out "$OUTPUT_DIR/data_final"

# Deleting unnecessary files
find $OUTPUT_DIR -type f -name "*data_sex_check*" -delete

# Calculating PCA
plink \
    --bfile "$OUTPUT_DIR/data_final" \
    --pca 5 \
    --out "$OUTPUT_DIR/pca5"
