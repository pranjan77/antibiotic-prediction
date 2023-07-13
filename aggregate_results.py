#!/usr/bin/env bash

import argparse
import pathlib
from typing import List

import pandas as pd


def update_results(results: pd.DataFrame):
    # Step1: Fitler results
    results_filtered = results[~results.function.isin({"antitumor", "antieuk"})]
    # Step 2: Update dataframe
    data_dict = dict()
    regions = set()
    for _, row in results_filtered.iterrows():
        if row.region not in regions:
            regions.add(row.region)
            data_item = dict()
            data_item["Genome"] = row.genome
            data_item["Region"] = row.region
            data_item[row.function] = row.probability
            data_dict[row.region] = data_item
        else:
            data_dict[row.region][row.function] = row.probability
    columns = [
        "Genome",
        "Region",
        "antibacterial",
        "antifungal",
        "antigrampos",
        "antigramneg",
    ]
    results_updated = pd.DataFrame(data_dict.values())[columns]
    return results_updated


def main(genome_paths: List[pathlib.Path], output_dir: pathlib.Path):
    partial_results = []
    genomes = [genome_path.stem for genome_path in genome_paths]
    for outputs in output_dir.iterdir():
        if outputs.is_dir() and outputs.stem in genomes:
            print("Processing results from " + str(outputs))
            partial_results_file = outputs / "prediction_results.csv"
            partial_results.append(pd.read_csv(partial_results_file))
    results = pd.concat(partial_results)
    results_updated = update_results(results)
    results_file_csv = output_dir / "aggregated_results.csv"
    results_file_html = output_dir / "aggregated_results.html"
    results_updated.to_csv(results_file_csv, sep=",", index=False)
    results_updated.to_html(results_file_html, index=False)


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument(
        "genomes",
        help="File paths of the genomes used for BGC function prediction",
        nargs="+",
    )
    PARSER.add_argument(
        "--output_dir", help="Directory for the output files", default="outputs"
    )
    ARGS = PARSER.parse_args()

    genome_paths = [pathlib.Path(genome) for genome in ARGS.genomes]
    output_dir = pathlib.Path(ARGS.output_dir)
    main(genome_paths, output_dir)
