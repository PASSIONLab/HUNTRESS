

# draw_tree() : visualize phylogenetic tree from a corresponding matrix
# Written by: Farid Rashidi Mehrabadi (frashidi@iu.edu)

import os
import time
import copy
import numpy as np
import pandas as pd
import networkx as nx
from datetime import datetime
import argparse
from argparse import ArgumentParser
import pygraphviz as pyg



def draw_tree(filename):
    add_cells = False
    change_edges_to_number = False
    combine = False

    from collections import Counter
    import pygraphviz as pyg
    import networkx as nx
    from networkx.drawing.nx_agraph import graphviz_layout, to_agraph

    def contains(col1, col2):
        for i in range(len(col1)):
            if not col1[i] >= col2[i]:
                return False
        return True

    df = pd.read_csv(filename, sep="\t", index_col=0)
    splitter_mut = "\n"
    matrix = df.values
    names_mut = list(df.columns)

    i = 0
    while i < matrix.shape[1]:
        j = i + 1
        while j < matrix.shape[1]:
            if np.array_equal(matrix[:, i], matrix[:, j]):
                matrix = np.delete(matrix, j, 1)
                x = names_mut.pop(j)
                names_mut[i] += splitter_mut + x
                j -= 1
            j += 1
        i += 1

    rows = matrix.shape[0]
    cols = matrix.shape[1]
    dimensions = np.sum(matrix, axis=0)
    indices = np.argsort(dimensions)
    dimensions = np.sort(dimensions)
    names_mut = [names_mut[indices[i]] for i in range(cols)]

    G = nx.DiGraph(dpi=300)
    G.add_node(cols)
    G.add_node(cols - 1)
    G.add_edge(cols, cols - 1, label=names_mut[cols - 1])
    node_mud = {}
    node_mud[names_mut[cols - 1]] = cols - 1

    i = cols - 2
    while i >= 0:
        if dimensions[i] == 0:
            break
        attached = False
        for j in range(i + 1, cols):
            if contains(matrix[:, indices[j]], matrix[:, indices[i]]):
                G.add_node(i)
                G.add_edge(node_mud[names_mut[j]], i, label=names_mut[i])
                node_mud[names_mut[i]] = i
                attached = True
                break
        if not attached:
            G.add_node(i)
            G.add_edge(cols, i, label=names_mut[i])
            node_mud[names_mut[i]] = i
        i -= 1

    clusters = {}
    for node in G:
        if node == cols:
            if not add_cells:
                G._node[node]["label"] = "<<b>zygote</b>>"
            G._node[node]["fontname"] = "Helvetica"
            G._node[node]["style"] = "rounded"
            G._node[node]["shape"] = "box"
            G._node[node]["margin"] = 0.05
            G._node[node]["pad"] = 0
            G._node[node]["width"] = 0
            G._node[node]["height"] = 0
            continue
        untilnow_mut = []
        sp = nx.shortest_path(G, cols, node)
        for i in range(len(sp) - 1):
            untilnow_mut += G.get_edge_data(sp[i], sp[i + 1])["label"].split(splitter_mut)
        untilnow_cell = df.loc[
            (df[untilnow_mut] == 1).all(axis=1)
            & (df[[x for x in df.columns if x not in untilnow_mut]] == 0).all(axis=1)
        ].index
        if len(untilnow_cell) > 0:
            clusters[node] = f'<b>{", ".join(untilnow_cell)}</b>'
        else:
            clusters[node] = "––"

        if add_cells:
            if "––" not in clusters[node]:
                G._node[node]["color"] = "#0000FF"
            else:
                G._node[node]["color"] = "gray70"
                G._node[node]["fontcolor"] = "gray70"
            G._node[node]["label"] = clusters[node]
        else:
            G._node[node]["label"] = ""
        G._node[node]["shape"] = "box"
        G._node[node]["fontname"] = "Helvetica"
        G._node[node]["style"] = "rounded"
        G._node[node]["margin"] = 0.05
        G._node[node]["pad"] = 0
        G._node[node]["width"] = 0
        G._node[node]["height"] = 0

    i = 1
    for k, v in clusters.items():
        if v == "––":
            clusters[k] = i * "––"
            i += 1

    outputpath = filename[: -len(".CFMatrix")]
    if add_cells:
        for node in G:
            if node != cols:
                num = 0
                paths = nx.shortest_path(G, source=cols, target=node)
                for i in range(len(paths) - 1):
                    x = paths[i]
                    y = paths[i + 1]
                    num += len(G[x][y]["label"].split(splitter_mut))
                G._node[node]["label"] = f"<[{node}]  " + G._node[node]["label"] + f"  ({num})>"
            else:
                G._node[node]["label"] = f"<[{node}]  <b>zygote</b>  (0)>"

    if change_edges_to_number:
        data = G.edges.data("label")
        # with open(f"{outputpath}.mutsAtEdges", "w") as fout:
        for u, v, l in data:
            ll = l.split(splitter_mut)
            # fout.write(f"[{u}]->[{v}]: {' '.join(ll)}\n")
            if "––" in G._node[v]["label"]:
                G.add_edge(u, v, label=f"  {len(ll)}  ", color="gray70", fontcolor="gray70")
            else:
                G.add_edge(u, v, label=f"  {len(ll)}  ")


    if combine:
        H = G.copy()

        for _ in range(len(G.nodes)):
            d_in = H.in_degree(H)
            d_out = H.out_degree(H)
            for node in H.nodes():
                if d_out[node] == 1 and d_in[node] == 1:
                    parent = [x for x in H.predecessors(node)][0]
                    child = [x for x in H.successors(node)][0]
                    if d_out[parent] < 2 and d_in[parent] == 1:
                        new_node = f'{parent}{node}'
                        new_edge = f"{int(H[parent][node]['label']) + int(H[node][child]['label'])}"

                        H = nx.contracted_nodes(H, parent, node, self_loops=False)
                        mapping = {parent: new_node}
                        H = nx.relabel_nodes(H, mapping)
                        H[new_node][child]['label'] = new_edge
                        break

        d_in = H.in_degree(H)
        d_out = H.out_degree(H)
        nodes = []
        for node in H.nodes():
            if d_out[node] == 0:
                nodes.append(node)

        for node in nodes:
            parent = [x for x in H.predecessors(node)][0]
            if d_out[parent] == 1 and d_in[parent] == 1:
                grandparent = [x for x in H.predecessors(parent)][0]

                new_node = f'{parent}{node}'
                new_edge = f"{int(H[grandparent][parent]['label']) + int(H[parent][node]['label'])}"

                H = nx.contracted_nodes(H, parent, node, self_loops=False)
                mapping = {parent: new_node}
                H = nx.relabel_nodes(H, mapping)
                H[grandparent][new_node]['label'] = new_edge

        mygraph = nx.drawing.nx_agraph.to_agraph(H)
    else:
        mygraph = nx.drawing.nx_agraph.to_agraph(G)

    mygraph.layout(prog="dot")
    mygraph.draw(f"{outputpath}.png")
