import networkx as nx
import random
import numpy as np
import pickle
import matplotlib.pyplot as plt
from networkx.utils import py_random_state
from networkx import empty_graph
import os
import multiprocessing
from multiprocessing import Process, Manager
import functools
import time
import math
import scipy.io as scio
from networkx.generators.classic import complete_graph, empty_graph, path_graph, star_graph

from utils_network import plot_temp_snapshot, graph_node_idx_permute, active_edge_generate_static

if __name__ == "__main__":
    path = "./results/"
    nodesnum = 100
    k = 6
    idx = 0
    snap_length = 100

    # %% ------- static network generate -------
    # graph = nx.random_graphs.barabasi_albert_graph(nodesnum, k // 2)
    # graph = nx.random_graphs.random_regular_graph(k, nodesnum)
    # mAdj = nx.to_numpy_array(graph)

    # np.save(path+'ba_static' + '_n' + str(nodesnum) + '_k' + str(k) + '.npy', mAdj)
    # np.save(path+'rr_static' + '_n' + str(nodesnum) + '_k' + str(k) + '.npy', mAdj)

    # %% ------- temporal network generate: randomly edge sample -------

    mAdj = np.load(path + 'ba_static' + '_n' + str(nodesnum) + '_k' + str(k) + '.npy')
    for p in np.round(np.arange(0.1, 1.0, 0.1), decimals=1):
        snap_mat = active_edge_generate_static(mAdj, p=p, snapshot=snap_length)
        file = path + 'ba_temp_snap' + str(snap_length) + '_n' + str(nodesnum) + '_k' + str(k) + "_p" + str(
            p) + '.npy'
        np.save(file, snap_mat)
