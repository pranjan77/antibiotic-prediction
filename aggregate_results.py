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
    genome_regions = set()
    for _, row in results_filtered.iterrows():
        genome_region = row.genome + "_" + row.region
        if genome_region not in genome_regions:
            genome_regions.add(genome_region)
            data_item = dict()
            data_item["Genome"] = row.genome
            data_item["Region"] = row.region
            data_item["Classifier"] = row.classifier
            data_item[row.function] = row.probability
            data_dict[genome_region] = data_item
        else:
            data_dict[genome_region][row.function] = row.probability
    columns = [
        "Genome",
        "Region",
        "Classifier",
        "antibacterial",
        "antifungal",
        "antigrampos",
        "antigramneg",
    ]
    results_updated = pd.DataFrame(data_dict.values())[columns]
    return results_updated


def generate_html_table(df: pd.DataFrame):
    """Display a pandas.DataFrame as jQuery DataTables"""
    import uuid

    # Generate random container name
    id_container = uuid.uuid1()
    output = """
<div id="datatable-container-{id_container}">
  <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.7.0/jquery.min.js"></script>
  <script type="text/javascript" src="https://cdn.datatables.net/1.13.5/js/jquery.dataTables.min.js"></script>
  <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.13.5/css/jquery.dataTables.min.css"/>
  <script type="text/javascript">
    $(document).ready( function () {{
        $('#BGCtable').DataTable();
    }});
  </script>
  <!-- Insert table below -->
  {table}
</div>
    """.format(
        id_container=id_container,
        table=df.to_html(
            index=False,
            table_id="BGCtable",
            classes="display",
            float_format=lambda x: "{:.3f}".format(x),
        ),
    )
    return output


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
    results_updated.to_csv(
        results_file_csv,
        sep=",",
        index=False,
        float_format=lambda x: "{:.3f}".format(x),
    )
    with open(results_file_html, "w") as fid:
        results_html = generate_html_table(results_updated)
        fid.write(results_html)


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
