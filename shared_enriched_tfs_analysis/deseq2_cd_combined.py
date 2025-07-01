# finding shared enriched TFs between DESeq2 and CD methods

import pandas as pd
import numpy as np
import json
import urllib.parse
from urllib.parse import quote
import requests
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from deseq2_analysis.DESeq2_TF_analysis import up_gene_list, down_gene_list, run_chea_kg, top_tfs

from characteristic_direction_analysis.CD_TF_analysis import cd_up_genes, cd_down_genes, run_chea_kg, top_tfs
from characteristic_direction_analysis.CD_TF_connections import fetch_chea_kg_data, get_tf_node_info, get_tf_edge_info

def shared_tfs(deseq2_input_csvs, cd_input_csvs, num_tfs):
    shared_tf_time_dict = {}
    num_comparisons = len(deseq2_input_csvs)
    for i in range(num_comparisons):
        deseq2_up_tfs = top_tfs(up_gene_list(deseq2_input_csvs[i]), num_tfs)
        deseq2_down_tfs = top_tfs(down_gene_list(deseq2_input_csvs[i]), num_tfs)
        cd_up_tfs = top_tfs(cd_up_genes(cd_input_csvs[i]), num_tfs)
        cd_down_tfs = top_tfs(cd_down_genes(cd_input_csvs[i]), num_tfs)

        shared_up_tfs = list(set(deseq2_up_tfs).intersection(set(cd_up_tfs)))
        shared_down_tfs = list(set(deseq2_down_tfs).intersection(set(cd_down_tfs)))
        shared_tf_time_dict[i] = (shared_up_tfs, shared_down_tfs)
    return shared_tf_time_dict


def create_tf_time_series_graph(shared_tf_time_dict, num_comparisons, filename):
    """
    Creates network of TFs given differentially expressed genes at each time point.
    """
    subnetwork = {"nodes": [], "edges": []}
    for time in shared_tf_time_dict.keys():
        up_tfs = shared_tf_time_dict[time][0]
        down_tfs = shared_tf_time_dict[time][1]

        for tf in up_tfs:
            node_info = get_tf_node_info(tf)
            node_info["data"]["color"] = "#80eaff"
            node_info["data"]["id"] = f"{node_info["data"]["label"]}_up_{time}"
            # node_info["data"]["label"] = f"{node_info["data"]["label"]}_up_{time}"
            subnetwork["nodes"].append(node_info)

        for tf in down_tfs:
            node_info = get_tf_node_info(tf)
            node_info["data"]["color"] = "#ff8a80"
            node_info["data"]["id"] = f"{node_info["data"]["label"]}_down_{time}"
            # node_info["data"]["label"] = f"{node_info["data"]["label"]}_down_{time}"
            subnetwork["nodes"].append(node_info)
    # print(subnetwork)
    # print(tf_time_dict)

    for i in range(num_comparisons-1):
        for j, source_tf_list in enumerate(shared_tf_time_dict[i]):
            for k, target_tf_list in enumerate(shared_tf_time_dict[i+1]):
                for source_tf in source_tf_list:
                    for target_tf in target_tf_list:
                        if source_tf != target_tf:
                            data = fetch_chea_kg_data(source_tf, target_tf)
                            if len(data["nodes"]) == 2 and len(data["edges"]) == 1:
                                if data["edges"][0]["data"]["source_label"] == source_tf and \
                                    data["edges"][0]["data"]["target_label"] == target_tf:
                                    if j == 0 and k == 0:
                                        data["edges"][0]["data"]["source"] = f"{source_tf}_up_{i}"
                                        # data["edges"][0]["data"]["source_label"] = f"{source_tf}_up_{i}"
                                        data["edges"][0]["data"]["target"] = f"{target_tf}_up_{i+1}"
                                        # data["edges"][0]["data"]["target_label"] = f"{target_tf}_up_{i+1}"
                                    elif j == 0 and k == 1:
                                        data["edges"][0]["data"]["source"] = f"{source_tf}_up_{i}"
                                        # data["edges"][0]["data"]["source_label"] = f"{source_tf}_up_{i}"
                                        data["edges"][0]["data"]["target"] = f"{target_tf}_down_{i+1}"
                                        # data["edges"][0]["data"]["target_label"] = f"{target_tf}_down_{i+1}"
                                    elif j == 1 and k == 0:
                                        data["edges"][0]["data"]["source"] = f"{source_tf}_down_{i}"
                                        # data["edges"][0]["data"]["source_label"] = f"{source_tf}_down_{i}"
                                        data["edges"][0]["data"]["target"] = f"{target_tf}_up_{i+1}"
                                        # data["edges"][0]["data"]["target_label"] = f"{target_tf}_up_{i+1}"
                                    elif j == 1 and k == 1:
                                        data["edges"][0]["data"]["source"] = f"{source_tf}_down_{i}"
                                        # data["edges"][0]["data"]["source_label"] = f"{source_tf}_down_{i}"
                                        data["edges"][0]["data"]["target"] = f"{target_tf}_down_{i+1}"
                                        # data["edges"][0]["data"]["target_label"] = f"{target_tf}_down_{i+1}"
                                    subnetwork["edges"].append(data["edges"][0])
                                    print("edge added, case 1")
                        else:
                            if len(source_tf_list) != 0:
                                try:
                                    edge = get_tf_edge_info(source_tf)
                                    if j == 0 and k == 0:
                                        edge["data"]["source"] = f"{source_tf}_up_{i}"
                                        # edge["data"]["source_label"] = f"{source_tf}_up_{i}"
                                        edge["data"]["target"] = f"{target_tf}_up_{i+1}"
                                        # edge["data"]["target_label"] = f"{target_tf}_up_{i+1}"
                                    elif j == 0 and k == 1:
                                        edge["data"]["source"] = f"{source_tf}_up_{i}"
                                        # edge["data"]["source_label"] = f"{source_tf}_up_{i}"
                                        edge["data"]["target"] = f"{target_tf}_down_{i+1}"
                                        # edge["data"]["target_label"] = f"{target_tf}_down_{i+1}"
                                    elif j == 1 and k == 0:
                                        edge["data"]["source"] = f"{source_tf}_down_{i}"
                                        # edge["data"]["source_label"] = f"{source_tf}_down_{i}"
                                        edge["data"]["target"] = f"{target_tf}_up_{i+1}"
                                        # edge["data"]["target_label"] = f"{target_tf}_up_{i+1}"
                                    elif j == 1 and k == 1:
                                        edge["data"]["source"] = f"{source_tf}_down_{i}"
                                        # edge["data"]["source_label"] = f"{source_tf}_down_{i}"
                                        edge["data"]["target"] = f"{target_tf}_down_{i+1}"
                                        # edge["data"]["target_label"] = f"{target_tf}_down_{i+1}"
                                    subnetwork["edges"].append(edge)
                                    print("edge added, case 2")
                                except ValueError:
                                    continue
    with open(f"{filename}.json", "w") as outfile:
        json.dump(subnetwork, outfile, indent=4)
    print("FINISHED")

# deseq2_input_csvs = ["deseq2_degs_1v0.csv", "deseq2_degs_3v1.csv", "deseq2_degs_6v3.csv", "deseq2_degs_12v6.csv", "deseq2_degs_24v12.csv"]
# cd_input_csvs = ["cd_degs_1v0.csv", "cd_degs_3v1.csv", "cd_degs_6v3.csv", "cd_degs_12v6.csv", "cd_degs_24v12.csv"]
# shared_tf_time_dict = shared_tfs(deseq2_input_csvs, cd_input_csvs, 10)
# print(shared_tf_time_dict)
# print(create_tf_time_series_graph(shared_tf_time_dict, 5, "shared_top_10_tfs_adjacent_time_pts"))

deseq2_input_csvs = ["deseq2_degs_1v0.csv", "deseq2_degs_3v0.csv", "deseq2_degs_6v0.csv", "deseq2_degs_12v0.csv", "deseq2_degs_24v0.csv"]
cd_input_csvs = ["cd_degs_1v0.csv", "cd_degs_3v0.csv", "cd_degs_6v0.csv", "cd_degs_12v0.csv", "cd_degs_24v0.csv"]
shared_tf_time_dict = shared_tfs(deseq2_input_csvs, cd_input_csvs, 10)
print(shared_tf_time_dict)
# print(create_tf_time_series_graph(shared_tf_time_dict, 5, "shared_top_10_tfs_compare_w_time_pt_0"))
