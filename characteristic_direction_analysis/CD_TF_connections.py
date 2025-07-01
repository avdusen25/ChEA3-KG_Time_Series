# determining the connections between TFs

import pandas as pd
import numpy as np
import json
import urllib.parse
from urllib.parse import quote
import requests
from .CD_TF_analysis import cd_up_genes
from .CD_TF_analysis import cd_down_genes
from .CD_TF_analysis import run_chea_kg
from .CD_TF_analysis import top_tfs

def fetch_chea_kg_data(start_tf, end_tf):
    """
    Outputs JSON data for shortest path connecting two TFs.
    """
    base_url = "https://chea-kg.maayanlab.cloud/api/knowledge_graph"

    query_filter = {
        "start": "Transcription Factor",
        "start_field": "label",
        "start_term": start_tf,
        "end": "Transcription Factor",
        "end_field": "label",
        "end_term": end_tf
    }

    encoded_filter = urllib.parse.quote(str(query_filter).replace("'", '"'))
    full_url = f"{base_url}?filter={encoded_filter}"

    response = requests.get(full_url)
    response.raise_for_status()
    return response.json()


def get_tf_node_info(tf_label):
    """
    Gets node information associated with single TF.
    """
    base_url = "https://chea-kg.maayanlab.cloud/api/knowledge_graph"

    query_filter = {
        "start": "Transcription Factor",
        "start_field": "label",
        "start_term": tf_label
    }
    encoded_filter = urllib.parse.quote(str(query_filter).replace("'", '"'))
    full_url = f"{base_url}?filter={encoded_filter}"
    # print(full_url)

    response = requests.get(full_url)
    response.raise_for_status()
    data = response.json()

    for node in data["nodes"]:
        if node["data"]["label"] == tf_label:
            return node
    raise ValueError(f"Node for TF '{tf_label}' not found in response.")


def get_tf_edge_info(tf_label):
    """
    Returns edge info if TF is autoregulatory.
    """
    base_url = "https://chea-kg.maayanlab.cloud/api/knowledge_graph"

    query_filter = {
        "start": "Transcription Factor",
        "start_field": "label",
        "start_term": tf_label
    }
    encoded_filter = urllib.parse.quote(str(query_filter).replace("'", '"'))
    full_url = f"{base_url}?filter={encoded_filter}"

    response = requests.get(full_url)
    response.raise_for_status()
    data = response.json()

    for edge in data["edges"]:
        if edge["data"]["source_label"] == tf_label and edge["data"]["target_label"] == tf_label:
            return edge
    raise ValueError(f"Self-edge for TF '{tf_label}' not found in response.")


def create_tf_time_series_graph(input_degs, num_comparisons, filename):
    """
    Creates network of TFs given differentially expressed genes at each time point.
    """
    subnetwork = {"nodes": [], "edges": []}
    time = 0
    tf_time_dict = {}

    for file in input_degs:
        up_tfs = top_tfs(cd_up_genes(file))
        down_tfs = top_tfs(cd_down_genes(file))
        tf_time_dict[time] = (up_tfs, down_tfs)

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

        time += 1
    # print(subnetwork)
    # print(tf_time_dict)

    for i in range(num_comparisons-1):
        for j, source_tf_list in enumerate(tf_time_dict[i]):
            for k, target_tf_list in enumerate(tf_time_dict[i+1]):
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

if __name__ == "__main__":
    # print(get_tf_node_info("PITX3"))

    # input_degs = ["cd_degs_1v0.csv",
    #             "cd_degs_3v1.csv",
    #             "cd_degs_6v3.csv",
    #             "cd_degs_12v6.csv",
    #             "cd_degs_24v12.csv"]
    # create_tf_time_series_graph(input_degs, 5, "cd_adjacent_time_pts_top_5_tfs")

    input_degs = ["cd_degs_1v0.csv",
                  "cd_degs_3v0.csv",
                  "cd_degs_6v0.csv",
                  "cd_degs_12v0.csv",
                  "cd_degs_24v0.csv"]
    create_tf_time_series_graph(input_degs, 5, "cd_compare_w_time_pt_0_top_5_tfs")
