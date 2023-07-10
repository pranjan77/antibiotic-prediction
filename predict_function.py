#!/usr/bin/env python3

import argparse
import pathlib
import pickle
from typing import List

from Bio import SeqIO
import numpy as np
import pandas as pd

from predict_function_tools import read_training_data, read_SSN_features
import readInputFiles


def collect_input_files(antismash_dir: pathlib.Path, rgi_dir: pathlib.Path):
    input_files = []
    antismash_files = list(antismash_dir.glob("*region*.gbk"))
    for antismash_file in antismash_files:
        as_filename = antismash_file.stem
        rgi_file = rgi_dir / (as_filename + "/" + as_filename + ".txt")
        input_files.append((antismash_file, rgi_file))
    return input_files


def read_classifiers(data_dir: pathlib.Path, classifier_list: List[str]):
    classifiers = dict()
    pickle_files = list((data_dir / "classifiers").glob("*.pkl"))
    for classifier_name in classifier_list:
        for pickle_file in pickle_files:
            fname = pickle_file.stem
            key = fname.rsplit("_", 1)[0]
            if fname.startswith(classifier_name):
                with open(pickle_file, "rb") as fid:
                    classifiers[key] = pickle.load(fid)
    return classifiers


def read_antismash_bgc(file: pathlib.Path):
    try:
        with open(file, "rU") as fid:
            record = SeqIO.read(fid, "genbank")
    except:
        raise ValueError("error reading antismash output file")
    as_features = record.features
    return as_features


def read_rgi_bgc(file: pathlib.Path):
    try:
        rgi_infile = open(file, "r")
    except:
        raise ValueError("error reading rgi output file")
    return rgi_infile


def run_classifiers(classifiers, test_features, antismash_bgc_file):
    predictions = []
    bgc_name = antismash_bgc_file.stem
    for cl_name, cl in classifiers.items():
        classifier_type, prediction_type = cl_name.split("_")
        probability = cl.predict_proba(test_features)
        prediction = {
            "bgc_name": bgc_name,
            "classifier_type": classifier_type,
            "prediction_type": prediction_type,
            "probability": probability[0, 1],
        }
        predictions.append(prediction)
        print(prediction)
    return predictions


def predict_function(
    antismash_bgc_file: pathlib.Path,
    rgi_bgc_file: pathlib.Path,
    data_dir: pathlib.Path,
    training_features,
    classifiers,
    antismash_version: int,
    rgi_version: int,
):
    as_features = read_antismash_bgc(antismash_bgc_file)
    rgi_infile = read_rgi_bgc(rgi_bgc_file)
    test_SSN_feature_matrix = read_SSN_features(data_dir, antismash_bgc_file)
    data_path = str(data_dir) + "/"
    test_features = readInputFiles.readInputFiles(
        as_features,
        antismash_version,
        rgi_infile,
        rgi_version,
        training_features,
        data_path,
        test_SSN_feature_matrix,
    )
    prediction_results = run_classifiers(classifiers, test_features, antismash_bgc_file)
    rgi_infile.close()
    return prediction_results


def main(
    antismash_dir: pathlib.Path,
    rgi_dir: pathlib.Path,
    data_dir: pathlib.Path,
    output_dir: pathlib.Path,
    classifier_list: List[str],
    antismash_version: int,
    rgi_version: int,
) -> None:
    data = []
    # Step 1: Get antismash and rgi files
    input_files = collect_input_files(antismash_dir, rgi_dir)
    # Step 2: Read classifier pickle files
    classifiers = read_classifiers(data_dir, classifier_list)
    # Step 3: For each BGC predict function and store data_items
    training_features = read_training_data(data_dir, antismash_version, rgi_version)
    for antismash_bgc_file, rgi_bgc_file in input_files:
        data_items = predict_function(
            antismash_bgc_file,
            rgi_bgc_file,
            data_dir,
            training_features,
            classifiers,
            antismash_version,
            rgi_version,
        )
        data.extend(data_items)
    df = pd.DataFrame(data)
    results_file = output_dir / "prediction_results.csv"
    df.to_csv(results_file, sep=",", index=False)
    return None


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument(
        "antismash_dir", help="Directory containing the antismash output"
    )
    PARSER.add_argument("rgi_dir", help="Directory containing the rgi output")
    PARSER.add_argument(
        "--data_dir", help="Directory containing the data (classifiers)", default="data"
    )
    PARSER.add_argument(
        "--output_dir", help="Directory for the output files", default="outputs"
    )
    PARSER.add_argument(
        "--classifiers",
        nargs="+",
        help="The types of classifiers to use",
        default=["tree"],
    )
    PARSER.add_argument(
        "--antismash_version",
        help="version of antismash used to generate antismash input file",
        default="5",
    )
    PARSER.add_argument(
        "--rgi_version",
        help="version of rgi used to generate antismash input file",
        default="5",
    )
    ARGS = PARSER.parse_args()

    antismash_dir = pathlib.Path(ARGS.antismash_dir)
    rgi_dir = pathlib.Path(ARGS.rgi_dir)
    data_dir = pathlib.Path(ARGS.data_dir)
    output_dir = pathlib.Path(ARGS.output_dir)
    classifier_list = ARGS.classifiers
    antismash_version = int(ARGS.antismash_version)
    rgi_version = int(ARGS.rgi_version)
    main(
        antismash_dir,
        rgi_dir,
        data_dir,
        output_dir,
        classifier_list,
        antismash_version,
        rgi_version,
    )
