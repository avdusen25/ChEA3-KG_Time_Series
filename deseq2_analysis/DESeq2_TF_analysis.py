# Analyzing DESeq2 DEGs from GSE271120 time series RNA-seq data

import pandas as pd
import numpy as np
import json
import urllib.parse
from urllib.parse import quote
import requests

def up_gene_list(input_csv, filename=None):
    """
    Outputs list of upregulated DEGs from DESeq2 results CSV.
    """
    df = pd.read_csv(input_csv)
    row_filter = df["log2FoldChange"] > 0
    filtered = df.loc[row_filter, df.columns[0]]
    gene_ids = list(filtered)

    up_list = []
    for gene in gene_ids:
        if "_" in gene:
            up_list.append(gene.split("_", 1)[1])
        else:
            up_list.append(gene)

    print(len(up_list))
    return up_list


def down_gene_list(input_csv, filename=None):
    """
    Outputs list of downregulated DEGs from DESeq2 results CSV.
    """
    df = pd.read_csv(input_csv)
    row_filter = df["log2FoldChange"] < 0
    filtered = df.loc[row_filter, df.columns[0]]
    gene_ids = list(filtered)

    down_list = []
    for gene in gene_ids:
        if "_" in gene:
            down_list.append(gene.split("_", 1)[1])
        else:
            down_list.append(gene)

    print(len(down_list))
    return down_list


def run_chea_kg(gene_list, num_tfs):
    """
    Outputs JSON of TF subnetwork best corresponding to input gene list.
    """
    CHEA_KG = 'https://chea-kg.maayanlab.cloud/api/enrichment'

    description = "Genes upregulated in SK-Mel-28 vs. primary melanocytes"
    payload = {
        'list': (None, "\n".join(gene_list)),
        'description': (None, description)
    }
    response=requests.post(f"{CHEA_KG}/addList", files=payload)
    data = json.loads(response.text)

    q = {
        'min_lib': 3, # minimum number of libraries that a TF must be ranked in
        'libraries': [
            {'library': "Integrated--meanRank", 'term_limit': num_tfs} # edit term_limit to change number of top-ranked TFs
        ],
        'limit':50, # controls number of edges returned - may cause issues with visualization if too large
        'userListId': data['userListId']
    }
    query_json=json.dumps(q)

    res = requests.post(CHEA_KG, data=query_json)
    if res.ok:
        data = json.loads(res.text)
        return data
    else:
        data = None
        return res.text


def top_tfs(gene_list, num_tfs=5):
    """
    Returns a list of the top N most enriched TFs corresponding to an input gene list.
    """
    enriched_tfs = run_chea_kg(gene_list, num_tfs)
    tfs_list = []
    for node in enriched_tfs["nodes"]:
        tfs_list.append(node["data"]["label"])
    return tfs_list


if __name__ == "__main__":
    print(top_tfs(up_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_3v0.csv"), 10))

    # print(up_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_1v0.csv", "degs_1v0"))
    # print(down_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_1v0.csv", "degs_1v0"))

    # print(up_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_3v1.csv", "degs_3v1"))
    # print(down_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_3v1.csv", "degs_3v1"))

    # print(up_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_6v3.csv", "degs_6v3"))
    # print(down_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_6v3.csv", "degs_6v3"))``

    # print(up_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_12v6.csv", "degs_12v6"))
    # print(down_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_12v6.csv", "degs_12v6"))

    # print(up_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_24v12.csv", "degs_24v12"))
    # print(down_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_24v12.csv", "degs_24v12"))

    # print(up_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_3v0.csv", "degs_3v0"))
    # print(down_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_3v0.csv", "degs_3v0"))

    # print(up_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_6v0.csv", "degs_6v0"))
    # print(down_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_6v0.csv", "degs_6v0"))

    # print(up_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_12v0.csv", "degs_12v0"))
    # print(down_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_12v0.csv", "degs_12v0"))

    # print(up_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_24v0.csv", "degs_24v0"))
    # print(down_gene_list("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_degs_24v0.csv", "degs_24v0"))
