# determining the connections between TFs

import pandas as pd
import numpy as np
import json
import urllib.parse
from urllib.parse import quote
import requests
from DESeq2_TF_analysis import up_gene_list
from DESeq2_TF_analysis import down_gene_list
from DESeq2_TF_analysis import run_chea_kg
from DESeq2_TF_analysis import top_tfs

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

# print(get_tf_node_info("TP53"))


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
        up_tfs = top_tfs(up_gene_list(file))
        down_tfs = top_tfs(down_gene_list(file))
        tf_time_dict[time] = (up_tfs, down_tfs)

        for tf in up_tfs:
            node_info = get_tf_node_info(tf)
            node_info["data"]["color"] = "#80eaff"
            node_info["data"]["id"] = f"{node_info['data']['label']}_up_{time}"
            # node_info["data"]["label"] = f"{node_info["data"]["label"]}_up_{time}"
            subnetwork["nodes"].append(node_info)

        for tf in down_tfs:
            node_info = get_tf_node_info(tf)
            node_info["data"]["color"] = "#ff8a80"
            node_info["data"]["id"] = f"{node_info['data']['label']}_down_{time}"
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
    # data1 = fetch_chea_kg_data("TP53", "NFKB1")
    # print(len(data1['nodes']), len(data1['edges']))
    # print(data1)

    # data2 = fetch_chea_kg_data("NFKB1", "TP53")
    # print(len(data2['nodes']), len(data2['edges']))
    # print(data2)
    # print(data1 == data2)

    input_degs = ["/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_1v0.csv",
                  "/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_3v1.csv",
                  "/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_6v3.csv",
                  "/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_12v6.csv",
                  "/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_24v12.csv"]
    # create_tf_time_series_graph(input_degs, 5, "adjacent_time_pts_top_5_tfs")

    # time = 0
    # tf_time_dict = {}
    # for file in input_degs:
    #     up_tfs = top_tfs(up_gene_list(file), 10)
    #     down_tfs = top_tfs(down_gene_list(file), 10)
    #     tf_time_dict[time] = (up_tfs, down_tfs)
    #     time += 1
    # print(tf_time_dict)

    # top 5 tfs
    tf_time_dict = {0: (['KLF6', 'FOSB', 'JUN', 'FOS', 'NR4A3'], ['NR4A1', 'BHLHE40', 'ZNF395', 'FOSB', 'ZNF324']), 1: (['BHLHE40', 'ATF3', 'FOSB', 'SNAI1', 'RELB'], ['BHLHE40', 'ZBED3', 'JUN', 'PPARG', 'NR2F2']), 2: (['TEAD1', 'SP3', 'ZBTB38', 'FOXM1', 'UBP1'], ['MAFF', 'HMGN3', 'ATF3', 'GTF3A', 'ZNF511']), 3: (['STAT2', 'SP100', 'BATF2', 'TRAFD1', 'IRF9'], ['GATAD2A', 'E2F1', 'ZBED4', 'ZNF598', 'SRCAP']), 4: (['TEAD1', 'STAT2', 'NFKB2', 'CREB3L2', 'STAT1'], ['ZNF239', 'ZNF146', 'PRMT3', 'FOSL1', 'CEBPZ'])}
    # top 10 tfs
    tf_time_dict = {0: (['ATF3', 'KLF6', 'FOSB', 'JUN', 'SNAI1', 'NFIL3', 'FOS', 'NR4A3', 'FOSL1', 'EGR2'], ['NR4A1', 'BHLHE40', 'ZNF395', 'ATF3', 'FOSB', 'ZNF740', 'JUN', 'ZNF594', 'PRDM8', 'ZNF324']), 1: (['BHLHE40', 'ATF3', 'MYC', 'JUNB', 'FOSB', 'FOSL2', 'NFKB1', 'SNAI1', 'RELB', 'FOSL1'], ['BHLHE40', 'ATF3', 'CREBL2', 'KLF7', 'ZBED3', 'FOXA1', 'JUN', 'ZNF608', 'PPARG', 'NR2F2']), 2: (['HMGA2', 'TEAD1', 'SP3', 'PRDM4', 'ZBED4', 'TCF20', 'ZBTB38', 'FOXM1', 'UBP1', 'GLYR1'], ['MAFF', 'HMGN3', 'ATF3', 'FOSB', 'ZNF581', 'EGR1', 'JUN', 'GTF3A', 'ZNF207', 'ZNF511']), 3: (['STAT2', 'CREB3', 'JUN', 'SP100', 'BATF2', 'PPARG', 'TRAFD1', 'IRF1', 'IRF9', 'STAT3'], ['GATAD2A', 'E2F1', 'MYC', 'FOXK2', 'ZBED4', 'ZNF598', 'TCF3', 'SRCAP', 'FOXM1', 'FOSL1']), 4: (['TEAD1', 'STAT2', 'HIF1A', 'NFKB2', 'CREB3L2', 'NFKB1', 'STAT1', 'MAFK', 'ZNF697', 'STAT3'], ['MYC', 'ZNF239', 'ZNF146', 'TFDP1', 'PRMT3', 'CEBPG', 'ETV4', 'FOSL1', 'HMGA1', 'CEBPZ'])}


    input_degs = ["/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_1v0.csv",
                  "/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_3v0.csv",
                  "/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_6v0.csv",
                  "/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_12v0.csv",
                  "/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_24v0.csv"]
    # create_tf_time_series_graph(input_degs, 5, "deseq2_compare_w_time_pt_0_top_5_tfs")

    time = 0
    tf_time_dict = {}
    for file in input_degs:
        up_tfs = top_tfs(up_gene_list(file), 10)
        down_tfs = top_tfs(down_gene_list(file), 10)
        tf_time_dict[time] = (up_tfs, down_tfs)
        time += 1
    print(tf_time_dict)

    # top 10 tfs
    tf_time_dict = {0: (['ATF3', 'KLF6', 'FOSB', 'JUN', 'SNAI1', 'NFIL3', 'FOS', 'NR4A3', 'FOSL1', 'EGR2'], ['NR4A1', 'BHLHE40', 'ZNF395', 'ATF3', 'FOSB', 'ZNF740', 'JUN', 'ZNF594', 'PRDM8', 'ZNF324']), 1: (['ATF3', 'MYC', 'JUNB', 'FOSB', 'FOSL2', 'JUN', 'SNAI1', 'NFIL3', 'NR4A3', 'FOSL1'], ['TEAD1', 'BHLHE40', 'CREBL2', 'ZNF436', 'TFCP2L1', 'TCF7L2', 'CREB3L2', 'JUN', 'KLF9', 'SOX13']), 2: (['ADNP2', 'PRDM4', 'FOXK2', 'ZBED4', 'HIVEP1', 'TCF20', 'BAZ2A', 'ZNF697', 'SRCAP', 'FOSL1'], ['CREB3', 'ELF3', 'CREB3L4', 'JUN', 'ZNF580', 'PPARG', 'NR1H3', 'KLF2', 'IRF9', 'ZNF524']), 3: (['HMGA2', 'ZNF267', 'HIVEP1', 'ZBTB11', 'NFKB2', 'MGA', 'ZNF134', 'NFKB1', 'RLF', 'RELB'], ['E2F1', 'E2F7', 'MBD3', 'THAP4', 'THAP7', 'ZNF837', 'ZNF580', 'TFDP1', 'ZNF511', 'HMGA1']), 4: (['TEAD1', 'ZNF407', 'HIVEP1', 'MGA', 'TCF20', 'SMAD3', 'RFX7', 'ASH1L', 'NFAT5', 'NCOA2'], ['E2F1', 'DRAP1', 'THAP4', 'THAP7', 'MYC', 'ZNF598', 'ZNF692', 'GTF3A', 'TFDP1', 'HMGA1'])}
