# converting CSVs of DEGs to GMT format

# given CSV of DEGs from DESeq2
import pandas as pd
from DESeq2_TF_analysis import up_gene_list, down_gene_list, run_chea_kg, top_tfs

def csv_to_gmt(input_csv_list, comparisons):
    gmt_dict = {}
    for i, file in enumerate(input_csv_list):
        up_genes = up_gene_list(file)
        print(len(up_genes))
        down_genes = down_gene_list(file)
        print(len(down_genes))
        gmt_dict[f"{comparisons[i]} up genes"] = up_genes
        gmt_dict[f"{comparisons[i]} down genes"] = down_genes

    with open("./deseq2_compare_w_time_pt_0_degs.gmt", "w") as file:
        for s,t in gmt_dict.items():
            file.write(str(s) + "\t\t" + "\t".join(t) + "\n")
    print("FINISHED")

# adjacent time point comparisons
input_csv_list = ["/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_1v0.csv",
                  "/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_3v1.csv",
                  "/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_6v3.csv",
                  "/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_12v6.csv",
                  "/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_24v12.csv"]
comparisons = ["Hour 1 vs Hour 0", "Hour 3 vs Hour 1", "Hour 6 vs Hour 3", "Hour 12 vs Hour 6", "Hour 24 vs Hour 12"]
# print(csv_to_gmt(input_csv_list, comparisons))

input_csv_list = ["/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_1v0.csv",
                  "/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_3v0.csv",
                  "/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_6v0.csv",
                  "/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_12v0.csv",
                  "/Users/andrewjenkinsvandusen/Downloads/mount_sinai_internship/deseq2_analysis/deseq2_degs_csvs/deseq2_degs_24v0.csv"]
comparisons = ["Hour 1 vs Hour 0", "Hour 3 vs Hour 0", "Hour 6 vs Hour 0", "Hour 12 vs Hour 0", "Hour 24 vs Hour 0"]
print(csv_to_gmt(input_csv_list, comparisons))


def gmt_to_tf_time_dict(gmt_file):
    """
    Converts GMT file containing DEGs to tf_time_dict.
    """
    with open(gmt_file, 'r') as f:
        lines = f.readlines()

    temp_dict = {}
    for line in lines:
        tokens = line.split("\t\t")
        term = tokens[0]
        genes = [x.split(',')[0].strip() for x in tokens[1].split('\t')]
        temp_dict[term] = top_tfs(genes)
        print("enriched TFs found")

    comparisons = list(temp_dict.keys())
    new_comparisons = []
    for item in comparisons:
        comp = item.rsplit(' ', 2)[0]
        if comp not in new_comparisons:
            new_comparisons.append(comp)
    print(new_comparisons)

    j = 0
    tf_time_dict = {}
    for i in range(len(comparisons) // 2):
        tf_time_dict[i] = (temp_dict[comparisons[j]], temp_dict[comparisons[j+1]])
        j += 2
    print(tf_time_dict)
    return tf_time_dict, new_comparisons

# gmt_to_tf_time_dict("deseq2_analysis/deseq2_adj_time_pts_degs.gmt")
