#!/usr/bin/env bash

OUTPUT_DIR=outputs
GENOMES=$@

# Step 1: Run BGC function prediction
for GENOME in $GENOMES; do
    echo "Running BGC function prediction on $GENOME"
    bash predict_function.sh $GENOME $OUTPUT_DIR
done

# Step 2: Aggregate results
echo "Aggregating results"
python aggregate_results.py $GENOMES --output_dir $OUTPUT_DIR
