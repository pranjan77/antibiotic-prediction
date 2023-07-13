#!/usr/bin/env bash

import argparse
import pathlib
from typing import List

import pandas as pd


def main(genome_paths: List[pathlib.Path], output_dir: pathlib.Path):
    partial_results = []
    genomes = [genome_path.stem for genome_path in genome_paths]
    for outputs in output_dir.iterdir():
        if outputs.is_dir() and outputs.stem in genomes:
            print("Processing results from " + str(outputs))
            partial_results_file = outputs / "prediction_results.csv"
            partial_results.append(pd.read_csv(partial_results_file))
    results = pd.concat(partial_results)
    results_file_csv = output_dir / "aggregated_results.csv"
    results_file_html = output_dir / "aggregated_results.html"
    results.to_csv(results_file_csv, sep=",", index=False)
    results.to_html(results_file_html, index=False)


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
