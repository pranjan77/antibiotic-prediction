#!/usr/bin/env bash

GENOME=$1
GENOME_NAME=$(basename "${GENOME%%.*}")
OUTPUT_DIR=$2/$GENOME_NAME
mkdir -p $OUTPUT_DIR
NO_SSN=$3

# Step 1: Run antismash
source deactivate &>/dev/null
source activate antismash
if [ ! -d "$OUTPUT_DIR/antismash" ]; then
    echo "Running antismash"
    mkdir -p $OUTPUT_DIR/antismash
    bash run_antismash7.sh $GENOME $OUTPUT_DIR/antismash
else
    echo "Antismash data for $GENOME already exists"
fi

# Step 2: For each BGC get fasta and run RGI
source deactivate &>/dev/null
source activate rgi5
if [ ! -d "$OUTPUT_DIR/rgi" ]; then
    echo "Running RGI"
    mkdir -p $OUTPUT_DIR/rgi
    for f in $OUTPUT_DIR/antismash/*region*.gbk; do
        BGC=$(basename "$f" .gbk)
        mkdir $OUTPUT_DIR/rgi/$BGC
        python gbk2fasta.py "$f" $OUTPUT_DIR/rgi/"$BGC".fna
        rgi main -i $OUTPUT_DIR/rgi/"$BGC".fna -o $OUTPUT_DIR/rgi/"$BGC"/"$BGC" --include_loose
    done
else
    echo "RGI data for $GENOME already exists"
fi

# Step 3: Run prediction script for each BGC
source deactivate &>/dev/null
source activate natural_product
if [ ! -f "$OUTPUT_DIR/prediction_results.csv" ]; then
    echo "Running BGC activity prediction"
    python predict_function.py $OUTPUT_DIR/antismash $OUTPUT_DIR/rgi \
        --data_dir data \
        --output_dir $OUTPUT_DIR \
        --no_SSN $NO_SSN \
        --classifiers tree \
        --antismash_version 5 \
        --rgi_version 5
else
    echo "BGC prediction data for $GENOME already exists"
fi

# TODO: Provide options for svm, log
