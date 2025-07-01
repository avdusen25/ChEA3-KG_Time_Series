# visualizing network of TFs given JSON data

from dash import Dash, html
import dash_cytoscape as cyto
import json

def visualize_network(data):
    stylesheet = [
                    {
                        "selector": 'node',
                        "style": {
                            'background-color': 'data(color)',
                            'border-color': 'data(borderColor)',
                            'border-width': 'data(borderWidth)',
                            'label': 'data(label)',
                            "text-valign": "center",
                            "text-halign": "center",
                            'width': 'mapData(node_type, 0, 1, 30, 80)',
                            'height': 'mapData(node_type, 0, 1, 30, 80)',
                        }
                    },
                    {
                        "selector": ".row-label",
                            "style": {
                            'background-color': '#ffffff',
                            'border-width': 0,
                            'label': 'data(label)',
                            'width': 1,
                            'height': 1,
                            'font-size': 20,
                            'text-valign': 'center',
                            'text-halign': 'center',
                            'color': '#000000',
                            'text-wrap': 'wrap',
                            'text-max-width': 100
                        }
                    },
                    {
                        "selector": 'edge',
                        "style": {
                            "curve-style": "bezier",
                            'line-color': 'data(lineColor)',
                            'width': '4',
                            'arrow-scale': 1.5,
                            "text-rotation": "autorotate",
                            "text-margin-x": 0,
                            "text-margin-y": 0,
                            'font-size': '12px',
                            'target-arrow-shape': 'data(directed)',
                            'target-endpoint': 'outside-to-node',
                            'source-endpoint': 'outside-to-node',
                            'target-arrow-color': 'data(lineColor)'
                        }
                    }
                ]
    app = Dash()

    app.layout = html.Div([

            cyto.Cytoscape(
            id='enrichment_results',
            layout={'name': 'grid', 'cols': 11},
            style={'width': '100%', 'height': '1000px'},
            elements=data,
            stylesheet=stylesheet
        )
    ])

    app.run()

"#ff8a80"
"#80eaff"

# deseq2 adjacent time pts
# with open("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_adjacent_time_pts_top_5_tfs.json", "r") as input:
#     data = json.load(input)
# comparisons = ["Hour 1 vs Hour 0", "Hour 3 vs Hour 1", "Hour 6 vs Hour 3", "Hour 12 vs Hour 6", "Hour 24 vs Hour 12"]
# cols = 11
# nodes = data["nodes"]
# for i, comparison in enumerate(comparisons):
#     row_index = i * cols
#     label_node = {
#         'data': {
#             'id': f'row_label_{i}',
#             'label': comparison,
#             'color': '#ffffff',
#             'borderColor': '#ffffff',
#             'borderWidth': 0
#         },
#         'classes': 'row-label'
#     }
#     nodes.insert(row_index, label_node)
# visualize_network(data)

# deseq2 time point 0 comparison
# with open("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/deseq2_compare_w_time_pt_0_top_5_tfs.json", "r") as input:
#     data = json.load(input)
# comparisons = ["Hour 1 vs Hour 0", "Hour 3 vs Hour 0", "Hour 6 vs Hour 0", "Hour 12 vs Hour 0", "Hour 24 vs Hour 0"]
# cols = 11
# nodes = data["nodes"]
# for i, comparison in enumerate(comparisons):
#     row_index = i * cols
#     label_node = {
#         'data': {
#             'id': f'row_label_{i}',
#             'label': comparison,
#             'color': '#ffffff',
#             'borderColor': '#ffffff',
#             'borderWidth': 0
#         },
#         'classes': 'row-label'
#     }
#     nodes.insert(row_index, label_node)
# visualize_network(data)


# cd adjacent time pts, top 5 tfs
# with open("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/cd_adjacent_time_pts_top_5_tfs.json", "r") as input:
#     data = json.load(input)
# comparisons = ["Hour 1 vs Hour 0", "Hour 3 vs Hour 1", "Hour 6 vs Hour 3", "Hour 12 vs Hour 6", "Hour 24 vs Hour 12"]
# cols = 11
# nodes = data["nodes"]
# for i, comparison in enumerate(comparisons):
#     row_index = i * cols
#     label_node = {
#         'data': {
#             'id': f'row_label_{i}',
#             'label': comparison,
#             'color': '#ffffff',
#             'borderColor': '#ffffff',
#             'borderWidth': 0
#         },
#         'classes': 'row-label'
#     }
#     nodes.insert(row_index, label_node)
# visualize_network(data)

# cd time point 0 comparison, top 10 tfs
# with open("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/cd_compare_w_time_pt_0_top_10_tfs.json", "r") as input:
#     data = json.load(input)
# visualize_network(data)

# cd time point 0 comparison, top 5 tfs
# with open("/Users/andrewjenkinsvandusen/Downloads/mount sinai internship/cd_compare_w_time_pt_0_top_5_tfs.json", "r") as input:
#     data = json.load(input)
# comparisons = ["Hour 1 vs Hour 0", "Hour 3 vs Hour 0", "Hour 6 vs Hour 0", "Hour 12 vs Hour 0", "Hour 24 vs Hour 0"]
# cols = 11
# nodes = data["nodes"]
# for i, comparison in enumerate(comparisons):
#     row_index = i * cols
#     label_node = {
#         'data': {
#             'id': f'row_label_{i}',
#             'label': comparison,
#             'color': '#ffffff',
#             'borderColor': '#ffffff',
#             'borderWidth': 0
#         },
#         'classes': 'row-label'
#     }
#     nodes.insert(row_index, label_node)
# visualize_network(data)


def visualize_network_2(data):
    row_col_dict = {}
    for node in data["nodes"]:
        if "position" in node:
            continue
        node_id = node["data"]["id"]
        row_index = int(node_id.split("_")[-1])
        col_index = row_col_dict.get(row_index, 1)
        row_col_dict[row_index] = col_index + 1 # so that col indices don't interfere next time you have the same row index
        node["position"] = {"x": col_index * 150, "y": row_index * 150}

    stylesheet = [
                    {
                        "selector": 'node',
                        "style": {
                            'background-color': 'data(color)',
                            'border-color': 'data(borderColor)',
                            'border-width': 'data(borderWidth)',
                            'label': 'data(label)',
                            "text-valign": "center",
                            "text-halign": "center",
                            'width': 'mapData(node_type, 0, 1, 30, 80)',
                            'height': 'mapData(node_type, 0, 1, 30, 80)',
                        }
                    },
                    {
                        "selector": '.row-label',
                            "style": {
                            'background-color': '#ffffff',
                            'border-width': 0,
                            'label': 'data(label)',
                            'width': 1,
                            'height': 1,
                            'font-size': 20,
                            'text-valign': 'center',
                            'text-halign': 'center',
                            'color': '#000000',
                            'text-wrap': 'wrap',
                            'text-max-width': 100
                        }
                    },
                    {
                        "selector": 'edge',
                        "style": {
                            "curve-style": "bezier",
                            'line-color': 'data(lineColor)',
                            'width': '4',
                            'arrow-scale': 1.5,
                            "text-rotation": "autorotate",
                            "text-margin-x": 0,
                            "text-margin-y": 0,
                            'font-size': '12px',
                            'target-arrow-shape': 'data(directed)',
                            'target-endpoint': 'outside-to-node',
                            'source-endpoint': 'outside-to-node',
                            'target-arrow-color': 'data(lineColor)'
                        }
                    }
                ]
    app = Dash()

    app.layout = html.Div([

            cyto.Cytoscape(
            id='enrichment_results',
            layout={'name': 'preset'},
            style={'width': '100%', 'height': '1000px'},
            elements=data,
            stylesheet=stylesheet
        )
    ])

    app.run()

# shared top 10 tfs, comparing adjacent time pts
# with open("shared_top_10_tfs_adjacent_time_pts.json", "r") as input:
#     data = json.load(input)
# comparisons = ["Hour 1 vs Hour 0", "Hour 3 vs Hour 1", "Hour 6 vs Hour 3", "Hour 12 vs Hour 6", "Hour 24 vs Hour 12"]
# nodes = data["nodes"]
# for i, comparison in enumerate(comparisons):
#     label_node = {
#         'data': {
#             'id': f'row_label_{i}',
#             'label': comparison,
#             'color': '#ffffff',
#             'borderColor': '#ffffff',
#             'borderWidth': 0
#         },
#         'classes': 'row-label',
#         'position': {"x": 0, "y": i * 150}
#     }
#     nodes.append(label_node)
# visualize_network_2(data)

# shared top 10 tfs, comparing w time pt 0
with open("shared_top_10_tfs_compare_w_time_pt_0.json", "r") as input:
    data = json.load(input)
comparisons = ["Hour 1 vs Hour 0", "Hour 3 vs Hour 0", "Hour 6 vs Hour 0", "Hour 12 vs Hour 0", "Hour 24 vs Hour 0"]
nodes = data["nodes"]
for i, comparison in enumerate(comparisons):
    label_node = {
        'data': {
            'id': f'row_label_{i}',
            'label': comparison,
            'color': '#ffffff',
            'borderColor': '#ffffff',
            'borderWidth': 0
        },
        'classes': 'row-label',
        'position': {"x": 0, "y": i * 150}
    }
    nodes.append(label_node)
visualize_network_2(data)
