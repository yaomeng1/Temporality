import networkx as nx
import random
import numpy as np
import pickle
import matplotlib.pyplot as plt
import scipy.io as scio
from networkx.utils import py_random_state
from networkx import empty_graph
import os
import multiprocessing
from multiprocessing import Process, Manager
import functools
import time
import math
from utils import temp_edge_dict, temp_nbr_dict, rand_pick, rand_pick_list


def mat_from_edge(edge_list_dict):
    snapnum = len(edge_list_dict)
    n = 200
    subnet_mat = np.zeros([snapnum, n, n], dtype=np.int)
    for snap in range(snapnum):
        edge_list = edge_list_dict[snap]
        for nodex, nodey in edge_list:
            subnet_mat[snap, nodex, nodey] = 1
            subnet_mat[snap, nodey, nodex] = 1

    return subnet_mat


def subnet_generate(static_matrix, snapshot=100, p=0.3):
    nodesnum = static_matrix.shape[0]
    graph = nx.from_numpy_matrix(static_matrix)
    edge_list = [edge for edge in graph.edges()]
    subnet_edge_list = {}
    sub_matrix_list = []
    for i in range(snapshot):
        sub_matrix = np.zeros([nodesnum, nodesnum], dtype=np.int)
        subnet_edge_list[i] = [edge for edge in edge_list if rand_pick(p)]
        # print(len(subnet_edge_list[i]))
        for nodex, nodey in subnet_edge_list[i]:
            sub_matrix[nodex, nodey] = 1
            sub_matrix[nodey, nodex] = 1
        sub_matrix_list.append(sub_matrix)
    return sub_matrix_list


def k_mean_var(mAdj):
    k_array = np.sum(mAdj, axis=1)
    return np.var(k_array)


def snap_struct(snap_mat):
    graph = nx.from_numpy_matrix(snap_mat)
    # comp_list = [c for c in nx.connected_components(graph)]
    comp_list = [c for c in sorted(nx.connected_components(graph), key=len, reverse=True)]

    return comp_list


def temp_struct(subnet_mat):
    """

    :param subnet_mat:
    :return: the max component scale; max component k mean; snapshot k mean;
    max component k variance
    """
    max_comp_list = []
    comp_k_list = []
    comp_kvar_list = []
    k_mean_list = []
    for snap in range(subnet_mat.shape[0]):
        snap_mat = subnet_mat[snap, :, :]
        k_mean = np.sum(snap_mat) / snap_mat.shape[0]
        k_mean_list.append(k_mean)
        comp_list = snap_struct(snap_mat)  # list of nodes set
        max_comp_set = comp_list[0]
        max_comp_list.append(len(max_comp_set))
        comp_k_array = np.sum(snap_mat[np.array(list(max_comp_set)), :], axis=1)

        comp_k_list.append(np.mean(comp_k_array))
        comp_kvar_list.append(np.var(comp_k_array))
    return np.array(max_comp_list), np.array(comp_k_list), \
           np.array(comp_kvar_list), np.array(k_mean_list)


def write_max_component(temp_mat):
    snapshot = temp_mat.shape[0]

    temp_mat_max_comp = np.zeros(temp_mat.shape, dtype=np.float_)
    for i in range(snapshot):
        graph = nx.from_numpy_array(temp_mat[i, :, :])
        comp_list = [c for c in sorted(nx.connected_components(graph), key=len, reverse=True)]
        max_comp_set = np.array(list(comp_list[0]))
        sub_mat_row = temp_mat[i, max_comp_set, :]
        temp_mat_max_comp[i, :len(max_comp_set), :len(max_comp_set)] = sub_mat_row[:, max_comp_set]

    return temp_mat_max_comp


def check_isolate(subnet):
    snapnum = subnet.shape[0]
    nodesnum = subnet.shape[1]
    iso_set = set(range(nodesnum))
    for i in range(snapnum):
        snap_mat = subnet[i, :, :]
        snap_iso_set = set(np.where(np.sum(snap_mat, axis=1) == 0)[0])
        iso_set = iso_set & snap_iso_set

    return iso_set


if __name__ == "__main__":

    # path = "E:\\pycharmproject\\work1\\TemporalEvolutionTheory\\"
    with open("./comp_fc_fix/sf_200_k6_4.pk", 'rb') as f:
        static_matrix = pickle.load(f)
    # with open("./comp_fc_fix/sf_200_k6_4.pk", 'rb') as f:
    #     static_matrix = pickle.load(f)
    for p in np.arange(0.5, 0.525, 0.05):
        sub_matrix_list = subnet_generate(static_matrix, snapshot=100, p=p)
        temp_mat = np.stack(sub_matrix_list)
        temp_mat_giant = write_max_component(temp_mat)
        temp_mat_giant = temp_mat_giant.astype(np.float_)

        scio.savemat("./sf_200_k6_4_p" + str(1000 * p)[0:3] + "_snapmatrix.mat", {"matrix_snap": temp_mat_giant})

