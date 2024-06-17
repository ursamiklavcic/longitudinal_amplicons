from __future__ import division
import os, sys
import random
import numpy as np
import matplotlib.pyplot as plt
import pickle

import scipy.stats as stats
from scipy.stats import gamma
import scipy.spatial as spatial

#from macroecotools import obs_pred_rsquare
from itertools import combinations


n_timepoints = 16
T_range = list(range(n_timepoints))





def calculate_dissimilarity(afd_temporal, T):

    # Eq. 3 in main text of Zaoli and Grilli

    diss_all = []
    # number of dissimilarities decreases with lag
    for t in range(len(afd_temporal) - T):

        n_t = afd_temporal[t]
        n_t_plus_T = afd_temporal[t + T]

        d_plus = n_t + n_t_plus_T
        d_minus = n_t - n_t_plus_T

        diss = ((d_minus**2) - d_plus)/(d_plus * (d_plus-1))
        diss_all.append(diss)

    mean_diss = np.mean(diss_all)

    return mean_diss



def calculate_dissimilarity_between_communities(afd_temporal_1, afd_temporal_2, T):

    # Form of Eq. 3 for between two ASVs within the same community.

    diss_all = []
    # number of dissimilarities decreases with lag
    for t in range(len(afd_temporal_1) - T):

        n_1_t = afd_temporal_1[t]
        n_2_t_plus_T = afd_temporal_2[t + T]

        d_plus = n_1_t + n_2_t_plus_T
        d_minus = n_1_t - n_2_t_plus_T

        diss = ((d_minus**2) - d_plus)/(d_plus * (d_plus-1))
        diss_all.append(diss)

    mean_diss = np.mean(diss_all)

    return mean_diss






# afd_temporal_all = list where each element is an array afd_temporal
# afd_temporal is an array of the read counts over time of a single ASV in a given community
# afd_temporal_all is of length S
# asv_names is of length S

# you'll need code to write this
afd_temporal_all = []
asv_name_all = []


diss_dict = {}
diss_dict['within_asv'] = {}
diss_dict['between_asv']  = {}
for afd_temporal_idx, afd_temporal in enumerate(afd_temporal_all):

    # remove if there's a zero
    afd_temporal = np.asarray(afd_temporal)
    # focus on ASVs that are present in all sampled
    # skip if there are any zeros (absences)
    if sum(afd_temporal==0) > 0:
        continue


    asv_name = asv_name_all[afd_temporal_idx]

    # Eq. 3 in main text of Zaoli and Grilli
    T_all = [calculate_dissimilarity(afd_temporal, t) for t in T_range]
    T_all = np.asarray(T_all)


    if asv_name not in diss_dict:
        diss_dict[asv_name] = {}

    diss_dict['within_asv'][asv_name]['afd_temporal'] = afd_temporal
    diss_dict['within_asv'][asv_name]['dissimilarity_within'] = T_all





# get pairs of ASVs
# go through species and get dissimilarity between replicate communities
asv_pair_all = list(combinations(list(diss_dict['within_asv'].keys()), 2))
diss_dict['between_asv']['asv_pair_all'] = asv_pair_all
diss_dict['between_asv']['asv_pair_all']['dissimilarity_between'] = []
for asv_pair in asv_pair_all:

    T_between_all = [calculate_dissimilarity_between_communities(diss_dict['within_asv'][asv_pair[0]]['afd_temporal'], diss_dict['within_asv'][asv_pair[1]]['afd_temporal'], t) for t in T_range]
    T_between_all = np.asarray(T_between_all)

    diss_dict['between_asv']['asv_pair_all']['dissimilarity_between'].append(T_between_all)



# plot the following 
# x-axis: T_range
# y-axis: mean of 'dissimilarity_within' over all ASVs for Fig. 2c
 


