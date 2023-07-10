#!/usr/bin/env bash

GENOME=$1
GENOME_FOLDER=$(basename $GENOME .fna)
OUTPUT_DIR=$2/$GENOME_FOLDER

# Step 1: Run antismash
bash run_antismash5.sh $GENOME $OUTPUT_DIR/antismash

# Step 2: For each BGC get fasta and run RGI
source deactivate
source activate rgi5
mkdir outputs/rgi
for f in $OUTPUT_DIR/antismash/*region*.gbk; do
    BGC=$(basename "$f" .gbk)
    mkdir $OUTPUT_DIR/rgi/$BGC
    python gbk2fasta.py "$f" $OUTPUT_DIR/rgi/"$BGC".fna
    rgi main -i $OUTPUT_DIR/rgi/"$BGC".fna -o $OUTPUT_DIR/rgi/"$BGC"/"$BGC" --include_loose
done

# Step 3: Run prediction script for each BGC
source deactivate
source activate natural_product
python predict_function.py $OUTPUT_DIR/antismash $OUTPUT_DIR/rgi \
    --data_dir data \
    --output_dir $OUTPUT_DIR \
    --classifiers tree \
    --antismash_version 5 \
    --rgi_version 5
