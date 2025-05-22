import networkx as nx
import random
import numpy as np
import pickle
import matplotlib.pyplot as plt
from numba import jit
import os
import multiprocessing
from multiprocessing import Process, Manager
import functools
import time
import math
from utils import temp_edge_dict, temp_nbr_dict, rand_pick, rand_pick_list


@jit(nopython=True)
def single_round(state_dict, payoff_dict, game_matrix, edge_mat):
    """
    play game on each group
    """

    for i in range(edge_mat.shape[0]):
        nodex = edge_mat[i, 0]
        nodey = edge_mat[i, 1]
        payoff_dict[nodex] += game_matrix[state_dict[nodex]][state_dict[nodey]]
        payoff_dict[nodey] += game_matrix[state_dict[nodey]][state_dict[nodex]]

    return payoff_dict


@jit(nopython=True)
def replicate_dynamic(state_dict, payoff_dict, nbr_mat, deg_array, nodesnum, w):
    """
    replicator dynamic after single round game: DB

    """

    update_node = np.random.choice(nodesnum)
    nbrs_num = deg_array[update_node]
    nbrs_array = nbr_mat[update_node][:nbrs_num]

    ## 如果选择到孤立节点，重新选择
    while nbrs_num == 0:
        update_node = np.random.choice(nodesnum)
        nbrs_num = deg_array[update_node]
        nbrs_array = nbr_mat[update_node][:nbrs_num]

    group_self = np.zeros(nbrs_num + 1, dtype=np.int_)
    group_self[:nbrs_num] = nbrs_array
    group_self[nbrs_num] = update_node

    fitness_group = 1 + w * payoff_dict[group_self]
    prob_array = fitness_group / np.sum(fitness_group)
    state_dict[update_node] = state_dict[rand_pick_list(group_self, prob_array)]

    return state_dict


def evolution(game_matrix, edge_mat_dict, nbr_mat_dict, nodesnum, w, snapshot_length, g):
    """
    whole process of evolution for 10000 times of generation
    """
    total_generation = int(1e8)

    payoff_dict = np.zeros(nodesnum)
    state_dict = np.zeros(nodesnum, dtype=np.int)
    coop_ini = random.choice(range(nodesnum))
    state_dict[coop_ini] = 1

    payoff_c_mean_array = np.zeros(total_generation, dtype=np.float_)
    payoff_mean_array = np.zeros(total_generation, dtype=np.float_)

    for time in range(total_generation):
        # if time % 1e4 == 0:
        #     print("time: ", time)
        idx = int(time / g) % snapshot_length
        payoff_dict = single_round(state_dict, payoff_dict, game_matrix, edge_mat_dict[idx])
        payoff_mean_array[time] = np.mean(payoff_dict)
        payoff_c_mean_array[time] = np.mean(payoff_dict[state_dict])
        nbr_mat, deg_array = nbr_mat_dict[idx]
        state_dict = replicate_dynamic(state_dict, payoff_dict, nbr_mat, deg_array, nodesnum, w)
        payoff_dict[:] = 0
        coord = np.sum(state_dict)
        if coord > nodesnum - 1:
            return 1
        if coord == 0:
            return 0

    return coord / nodesnum


def process(core, b, edge_mat_dict, nbr_mat_dict, nodesnum, g):
    w = 0.01
    snapshot_length = len(edge_mat_dict)
    game_matrix = np.zeros((2, 2))
    game_matrix[0][0] = 0  # P defect--defect
    game_matrix[0][1] = b  # T d-c
    game_matrix[1][0] = -1  # S
    game_matrix[1][1] = b - 1  # R

    repeat_time = int(1e4)
    repeat_array = np.zeros(repeat_time)

    for rep in range(repeat_time):
        coord_freq = evolution(game_matrix, edge_mat_dict, nbr_mat_dict,
                               nodesnum, w, snapshot_length, g)
        repeat_array[rep] = coord_freq

    return np.sum(repeat_array == 1) / (np.sum(repeat_array == 1) + np.sum(repeat_array == 0))


if __name__ == "__main__":
    with open("./rr_k6_p03_n200_edge_list_4.pk", 'rb') as f:
        edge_list_dict = pickle.load(f)
    # not the same thing: contact_matrix_edge is not symmetric, so edge on appear once in the matrix
    with open("./rr_k6_p03_n200_nbrs_dict_4.pk", 'rb') as f:
        nbrs_dict = pickle.load(f)

    edge_mat_dict = temp_edge_dict(edge_list_dict)
    nbr_mat_dict = temp_nbr_dict(nbrs_dict)
    nodesnum = 200
    g_list = [1e3, 2e3, 1e4]
    for g in g_list:
        fc_b_list = []
        print("g=", str(int(g)))
        if g == 1e3:
            b_para_list = np.around(np.arange(6.4, 7.4, 0.2), decimals=1)
        elif g == 2e3:
            b_para_list = np.around(np.arange(5.2, 6.2, 0.2), decimals=1)
        else:
            b_para_list = np.around(np.arange(4.2, 5.2, 0.2), decimals=1)
        for b_para in b_para_list:
            core_list = np.arange(20)  # 64-cpu core

            pool = multiprocessing.Pool()
            t1 = time.time()

            pt = functools.partial(process, b=b_para, edge_mat_dict=edge_mat_dict, nbr_mat_dict=nbr_mat_dict,
                                   nodesnum=nodesnum, g=g)
            coor_freq_list = pool.map(pt, core_list)

            coor_freq_core = sum(coor_freq_list) / len(coor_freq_list)

            pool.close()
            pool.join()
            t2 = time.time()
            print("Total time:" + (t2 - t1).__str__())
            print((b_para, coor_freq_core))
            fc_b_list.append(coor_freq_core)
        file = "./rr_200_k6_p03_im_g" + str(int(g)) + "_b" + str(min(b_para_list)) + "_" + str(
            max(b_para_list)) + "_2.pk"
        with open(file, 'wb') as f:
            pickle.dump(fc_b_list, f)
