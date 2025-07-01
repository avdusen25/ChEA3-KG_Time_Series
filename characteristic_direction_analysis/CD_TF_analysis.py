# finding DEGs in GSE271120 using characteristic direction

import numpy as np
import pandas as pd
# print(pd.__version__, np.__version__)
import json
from urllib.parse import quote
import requests
# from maayanlab_bioinformatics import normalization as norm
# from maayanlab_bioinformatics.dge import characteristic_direction
from dash import Dash, html
import dash_cytoscape as cyto

### Workflow given below for using CD to calculate DEGs (in actuality, ran it on a Colab notebook):
# preprocessed_data = pd.read_csv("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271120/suppl/GSE271120%5FRawCountFile%5Frsemgenes.CCBR1062.csv.gz")
# sample_columns = preprocessed_data.columns[1:]
# print(sample_columns)
# # # Index(['H00_S1', 'H00_S2', 'H00_S3', 'H01_S1', 'H01_S2', 'H03_S1', 'H03_S2',
# #        'H03_S3', 'H06_S1', 'H06_S2', 'H06_S3', 'H12_S1', 'H12_S2', 'H12_S3',
# #        'H24_S1', 'H24_S2', 'H24_S3'],
# #       dtype='object')
# preprocessed_data.set_index('gene_id', inplace=True)

# log_norm = norm.log10_normalize(preprocessed_data)
# q_norm = norm.quantile_normalize(log_norm) # quantile normalization
# norm_data = norm.zscore_normalize(q_norm.transpose()) # z-score normalization
# norm_data = norm_data.transpose()

# control_cols = ['H00_S1', 'H00_S2', 'H00_S3']
# case_cols = ['H01_S1', 'H01_S2']
# cases = preprocessed_data[case_cols]
# control = preprocessed_data[control_cols]
# all = pd.concat([cases, control], axis=1)
# cd = characteristic_direction(controls_mat=control, cases_mat=cases, calculate_sig=True)
# print(cd.head())


def cd_up_genes(input_csv):
    """
    Outputs list of upregulated genes determined by CD.
    """
    cd = pd.read_csv(input_csv, index_col=0)
    top_n = 600 # change this as desired
    highest_abs_expr = cd.loc[cd.abs().sort_values('CD-coefficient', ascending=False)[:top_n].index]
    up_genes = highest_abs_expr[highest_abs_expr['CD-coefficient'] > 0].dropna().index.to_list()

    up_list = []
    for gene in up_genes:
        gene_str = ""
        underscore = False
        for char in gene:
            if char == "_":
                underscore = True
            if underscore and char != "_":
                gene_str += char
        up_list.append(gene_str)

    print("Number of upregulated genes: ", len(up_genes))
    return up_list


def cd_down_genes(input_csv):
    """
    Outputs list of downregulated genes determined by CD.
    """
    cd = pd.read_csv(input_csv, index_col=0)
    top_n = 600 # change this as desired
    highest_abs_expr = cd.loc[cd.abs().sort_values('CD-coefficient', ascending=False)[:top_n].index]
    down_genes = highest_abs_expr[highest_abs_expr['CD-coefficient'] < 0].dropna().index.to_list()

    down_list = []
    for gene in down_genes:
        gene_str = ""
        underscore = False
        for char in gene:
            if char == "_":
                underscore = True
            if underscore and char != "_":
                gene_str += char
        down_list.append(gene_str)

    print("Number of downregulated genes: ", len(down_genes))
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
    # print(cd_up_genes("cd_degs_1v0.csv"))
    # print(cd_down_genes("cd_degs_1v0.csv"))

    print(top_tfs(cd_up_genes("cd_degs_24v0.csv"), 10))
    print(top_tfs(cd_down_genes("cd_degs_24v0.csv"), 10))


# groups
['H00_S1', 'H00_S2', 'H00_S3']
['H01_S1', 'H01_S2']
['H03_S1', 'H03_S2','H03_S3']
['H06_S1', 'H06_S2', 'H06_S3']
['H12_S1', 'H12_S2', 'H12_S3']
['H24_S1', 'H24_S2', 'H24_S3']
