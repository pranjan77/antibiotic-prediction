#!/usr/bin/env python3

import argparse
import pathlib
import pickle
from typing import List

from Bio import SeqIO
import numpy as np
import pandas as pd

import cluster_function_prediction_tools as tools
import readFeatureFiles
import readInputFiles
import SSN_tools


SSN_pfam_names = [
    "Thiolase, N-terminal domain",
    "ABC transporter",
    "Acyl transferase domain",
    "AAA domain",
    "ABC-2 family transporter protein",
    "Acyl-CoA dehydrogenase, C-terminal domain",
    "Acyl-CoA dehydrogenase, N-terminal domain",
    "Alcohol dehydrogenase GroES-like domain",
    "Alpha/beta hydrolase family",
    "Aminotransferase class I and II",
    "Beta-ketoacyl synthase, C-terminal domain",
    "Beta-ketoacyl synthase, N-terminal domain",
    "Cytochrome P450",
    "DegT/DnrJ/EryC1/StrS aminotransferase family",
    "Enoyl-(Acyl carrier protein) reductase",
    "Erythronolide synthase docking",
    "FAD binding domain",
    "Glycosyl transferase family 2",
    "Glycosyltransferase family 28 N-terminal domain",
    "Glycosyl transferases group 1",
    "Glycosyltransferase like family 2",
    "Glyoxalase/Bleomycin resistance protein/Dioxygenase superfamily",
    "KR domain",
    "Lanthionine synthetase C-like protein",
    "Major Facilitator Superfamily",
    "Methyltransferase small domain",
    "Methyltransferase domain",
    "NAD dependent epimerase/dehydratase family",
    "NDP-hexose 2,3-dehydratase",
    "O-methyltransferase",
    "Oxidoreductase family, C-terminal alpha/beta domain",
    "Oxidoreductase family, NAD-binding Rossmann fold",
    "Phosphopantetheine attachment site",
    "Polyketide cyclase / dehydrase and lipid transport",
    "Polyketide synthase dehydratase",
    "Protein of unknown function (DUF1205)",
    "short chain dehydrogenase",
    "SnoaL-like domain",
    "SpaB C-terminal domain",
    "Sugar (and other) transporter",
    "transcriptional_regulatory_protein,_c_terminal_domains",
    "Thioesterase superfamily",
    "ubiE/COQ5 methyltransferase family",
    "UDP-glucoronosyl and UDP-glucosyl transferase",
    "YcaO-like family",
    "Zinc-binding dehydrogenase",
    "pyridine_nucleotide-disulphide_oxidoreductase",
]


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


def read_training_data(
    data_dir: pathlib.Path, antismash_version: int, rgi_version: int
):
    data_path = str(data_dir) + "/"
    try:
        training_SSN_features = readFeatureFiles.readFeatureMatrix(
            data_path + "feature_matrices/SSN.csv"
        )
        if antismash_version == 4:
            training_pfam_features = readFeatureFiles.readFeatureMatrix(
                data_path + "feature_matrices/PFAM.csv"
            )
            training_smCOG_features = readFeatureFiles.readFeatureMatrix(
                data_path + "feature_matrices/SMCOG.csv"
            )
            # SSN_calc_features = readFeatureFiles.readFeatureMatrixFloat("gene_feature_matrices/test_compounds_SSN.csv")
            training_CDS_features = readFeatureFiles.readFeatureMatrix(
                data_path + "feature_matrices/CDS_motifs.csv"
            )

            training_pks_nrps_type_features = readFeatureFiles.readFeatureMatrix(
                data_path + "feature_matrices/pks_nrps_type.csv"
            )
            training_pk_signature_features = readFeatureFiles.readFeatureMatrix(
                data_path + "feature_matrices/pk_signature.csv"
            )
            training_pk_minowa_features = readFeatureFiles.readFeatureMatrix(
                data_path + "feature_matrices/pk_minowa.csv"
            )
            training_pk_consensus_features = readFeatureFiles.readFeatureMatrix(
                data_path + "feature_matrices/pk_consensus.csv"
            )

            training_nrp_stachelhaus_features = readFeatureFiles.readFeatureMatrix(
                data_path + "feature_matrices/nrp_stachelhaus.csv"
            )
            training_nrp_nrpspredictor_features = readFeatureFiles.readFeatureMatrix(
                data_path + "feature_matrices/nrp_nrpspredictor.csv"
            )
            training_nrp_pHMM_features = readFeatureFiles.readFeatureMatrix(
                data_path + "feature_matrices/nrp_pHMM.csv"
            )
            training_nrp_predicat_features = readFeatureFiles.readFeatureMatrix(
                data_path + "feature_matrices/nrp_predicat.csv"
            )
            training_nrp_sandpuma_features = readFeatureFiles.readFeatureMatrix(
                data_path + "feature_matrices/nrp_sandpuma.csv"
            )
        elif antismash_version == 5:
            training_pfam_features = readFeatureFiles.readFeatureMatrix(
                data_path + "feature_matrices/PFAM5.csv"
            )
            training_smCOG_features = readFeatureFiles.readFeatureMatrix(
                data_path + "feature_matrices/SMCOG5.csv"
            )
            training_CDS_features = readFeatureFiles.readFeatureMatrix(
                data_path + "feature_matrices/CDS_motifs5.csv"
            )
            training_pk_consensus_features = readFeatureFiles.readFeatureMatrix(
                data_path + "feature_matrices/pk_nrp_consensus5.csv"
            )

        if rgi_version == 3:
            training_card_features = readFeatureFiles.readFeatureMatrix(
                data_path + "feature_matrices/CARD_gene.csv"
            )
            used_resistance_genes_list = readFeatureFiles.readFeatureList(
                data_path + "feature_matrices/CARD_gene_list.txt"
            )
        elif rgi_version == 5:
            training_card_features = readFeatureFiles.readFeatureMatrix(
                data_path + "feature_matrices/CARD5_genes.csv"
            )
            used_resistance_genes_list = readFeatureFiles.readFeatureList(
                data_path + "feature_matrices/CARD5_gene_list.txt"
            )
    except:
        raise ValueError(
            "did not find file containing training data, please keep script located in directory downloaded from github"
        )

    training_features = np.concatenate(
        (training_pfam_features, training_card_features), axis=1
    )
    training_features = np.concatenate(
        (training_features, training_smCOG_features), axis=1
    )
    training_features = np.concatenate(
        (training_features, training_CDS_features), axis=1
    )
    training_features = np.concatenate(
        (training_features, training_SSN_features), axis=1
    )
    if antismash_version == 4:
        training_features = np.concatenate(
            (training_features, training_pks_nrps_type_features), axis=1
        )
        training_features = np.concatenate(
            (training_features, training_pk_signature_features), axis=1
        )
        training_features = np.concatenate(
            (training_features, training_pk_minowa_features), axis=1
        )
        training_features = np.concatenate(
            (training_features, training_pk_consensus_features), axis=1
        )
        training_features = np.concatenate(
            (training_features, training_nrp_stachelhaus_features), axis=1
        )
        training_features = np.concatenate(
            (training_features, training_nrp_nrpspredictor_features), axis=1
        )
        training_features = np.concatenate(
            (training_features, training_nrp_pHMM_features), axis=1
        )
        training_features = np.concatenate(
            (training_features, training_nrp_predicat_features), axis=1
        )
        training_features = np.concatenate(
            (training_features, training_nrp_sandpuma_features), axis=1
        )
    else:
        training_features = np.concatenate(
            (training_features, training_pk_consensus_features), axis=1
        )
    return training_features


def read_SSN_features(data_dir: pathlib.Path, antismash_bgc_file: pathlib.Path):
    data_path = str(data_dir) + "/"
    SSN_list = readFeatureFiles.readFeatureList(
        data_path + "feature_matrices/SSN_list.txt"
    )
    for i in range(0, len(SSN_list)):
        SSN_list[i] = SSN_list[i].replace("\r", "")
        SSN_list[i] = SSN_list[i].replace("\n", "")

    included_SSN_clusters = {}
    for pfam_name in SSN_list:
        base = pfam_name[0 : pfam_name.rfind("_")]
        if base not in included_SSN_clusters:
            included_SSN_clusters[base] = []
        numbering = pfam_name[pfam_name.rfind("_") + 1 : len(pfam_name)].replace(
            "\r", ""
        )
        included_SSN_clusters[base].append(numbering.replace("\n", ""))

    antismash_infilename = str(antismash_bgc_file)
    cluster_name = antismash_infilename
    blastp_path = "blastp"

    if "/" in cluster_name:
        cluster_name = cluster_name[cluster_name.rfind("/") + 1 : len(cluster_name)]
    cluster_name = cluster_name[0 : cluster_name.find(".gbk")]
    test_SSN_feature_matrix = SSN_tools.generateSSNFeatureMatrix(
        [antismash_infilename],
        SSN_pfam_names,
        SSN_list,
        included_SSN_clusters,
        blastp_path,
        cluster_name,
        data_path,
    )
    return test_SSN_feature_matrix


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
    classifiers,
    antismash_version: int,
    rgi_version: int,
):
    as_features = read_antismash_bgc(antismash_bgc_file)
    rgi_infile = read_rgi_bgc(rgi_bgc_file)
    training_features = read_training_data(data_dir, antismash_version, rgi_version)
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
    for antismash_bgc_file, rgi_bgc_file in input_files:
        data_items = predict_function(
            antismash_bgc_file,
            rgi_bgc_file,
            data_dir,
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
