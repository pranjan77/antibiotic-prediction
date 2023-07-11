#!/usr/bin/env bash

OUTPUT_DIR=outputs
NO_SSN="False"
GENOMES=$@

# Step 1: Run BGC function prediction
for GENOME in $GENOMES; do
    echo "Running BGC function prediction on $GENOME"
    bash predict_function.sh $GENOME $OUTPUT_DIR $NO_SSN
done

# Step 2: Aggregate results
echo "Aggregating results"
python aggregate_results.py $GENOMES --output_dir $OUTPUT_DIR
