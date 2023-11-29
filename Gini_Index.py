import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
from multiprocessing import Process
sys.path.append('/home/scientist617/Grad_Proj/Codes')
import re
from collections import defaultdict
from itertools import combinations
import networkx as nx
import community
from utils.cluster import cluster
from utils.recombinations import recombinations
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from pandas.io.formats.style import Styler
from utils.preprocessing import preprocessing
import os
import seaborn as sns
import igraph as ig
import scipy.stats as stats
import mplcursors
import plotly.express as px
import plotly.graph_objs as go
import plotly.io as pio
from numpy import linalg as LA
from Bio.PDB import PDBParser, Selection, NeighborSearch
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
from Bio.SeqUtils import seq3
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import nt_search
from Bio.Seq import MutableSeq
import logomaker
import statannot

cdrpos = {
    # http://www.bioinf.org.uk/abs/info.html#cdrid
    # light chain loop positions (chothia)
    'L': {
        'L1': ['26', '32'],
        'L2': ['50', '52'],
        'L3': ['91', '96'],
    },
    # heavy chain loop positions (chothia)
    'H': {
        'H1': ['26', '32'],
        'H2': ['52', '56'],
        'H3': ['96', '101']
    }
}
# %%
sjogren_path = '/home/scientist617/Grad_Proj/Data/new_sjogren_file/'
healthy_path = '/home/scientist617/Grad_Proj/Data/in_house_healthy/'



gene_list_all = ['IGHV1-18', 'IGHV1-24','IGHV1-69', 'IGHV3-23', 'IGHV3-43','IGHV3-53', 'IGHV3-72', 'IGHV4-34', 'IGHV5-51', 'IGHV6-1']

def gini_run_edges(gene):
    print(f"Running with input: {gene}")
    gini_sjogren = []
    gini_healthy = []
    for file in os.listdir(sjogren_path):
        df = preprocessing(pd.read_csv(os.path.join(sjogren_path, file), sep='\t'))()
        cl = cluster(df[df['v_call'] == gene])
        if len(df[df['v_call'] == gene]) < 10:
            continue
        else:
            _, temp = cl.clonal_expansion_index_by_edges()
            gini_sjogren.append(temp)
    for file in os.listdir(healthy_path):
        df = preprocessing(pd.read_csv(os.path.join(healthy_path, file), sep='\t'))()
        cl = cluster(df[df['v_call'] == gene])
        if len(df[df['v_call'] == gene]) < 10:
            continue
        else:
            _, temp = cl.clonal_expansion_index_by_edges()
            gini_healthy.append(temp)

    max_length = max(len(gini_sjogren), len(gini_healthy))

    if len(gini_sjogren)<max_length:
        gini_sjogren.extend([np.nan]*(max_length-len(gini_sjogren)))
    if len(gini_healthy)<max_length:
        gini_healthy.extend([np.nan]*(max_length-len(gini_healthy)))

    df_gini = pd.DataFrame({'Healthy': gini_healthy, 'Sjogren': gini_sjogren})
    plt.figure(figsize = (4, 4))
    sns.set(style="whitegrid")
    ax = sns.boxplot(data=df_gini, palette=['lightgrey', 'skyblue'], width=0.35)
    ax.set(ylim=(0, 1), xlim=(0, 1))
    sns.stripplot(data=df_gini, jitter=True, color=".3")

    statannot.add_stat_annotation(ax, data=df_gini, box_pairs=[('Healthy', 'Sjogren')], test='t-test_ind', text_format='star', verbose=2, loc="inside")

    # Customize the plot
    plt.ylabel('Gini Index (edges)')
    plt.title(f'{gene}')
    plt.tight_layout()
    plt.savefig(f'/home/scientist617/Grad_Proj/Data/Figures/Gini/Gini_index_nodes_{gene}.png', dpi=500)
    plt.close()
    print(f"file Gini_index_nodes_{gene}.png saved!")
def gini_run_duplicates(gene):
    print(f"Running with input: {gene}")
    gini_sjogren = []
    gini_healthy = []
    for file in os.listdir(sjogren_path):
        df = preprocessing(pd.read_csv(os.path.join(sjogren_path, file), sep='\t'))()
        cl = cluster(df[df['v_call'] == gene])
        if len(df[df['v_call'] == gene]) == 0:
            continue
        _, temp = cl.clonal_expansion_index_by_duplicates()
        gini_sjogren.append(temp)
    for file in os.listdir(healthy_path):
        df = preprocessing(pd.read_csv(os.path.join(healthy_path, file), sep='\t'))()
        cl = cluster(df[df['v_call'] == gene])
        if len(df[df['v_call'] == gene]) == 0:
            continue
        _, temp = cl.clonal_expansion_index_by_duplicates()
        gini_healthy.append(temp)

    max_length = max(len(gini_sjogren), len(gini_healthy))

    if len(gini_sjogren) < max_length:
        gini_sjogren.extend([np.nan] * (max_length - len(gini_sjogren)))
    if len(gini_healthy) < max_length:
        gini_healthy.extend([np.nan] * (max_length - len(gini_healthy)))

    df_gini = pd.DataFrame({'Healthy': gini_healthy, 'Sjogren': gini_sjogren})
    plt.figure(figsize=(4, 4))
    sns.set(style="whitegrid")
    ax = sns.boxplot(data=df_gini, palette=['lightgrey', 'skyblue'], width=0.35)
    ax.set(ylim=(0, 1), xlim=(0, 1))
    sns.stripplot(data=df_gini, jitter=True, color=".3")
    statannot.add_stat_annotation(ax, data=df_gini, box_pairs=[('Healthy', 'Sjogren')], test='t-test_ind',
                                  text_format='star', verbose=2, loc="inside")

    # Customize the plot
    plt.ylabel('Gini Index (duplicates)')
    plt.title(f'{gene}')
    plt.tight_layout()
    plt.savefig(f'/home/scientist617/Grad_Proj/Data/Figures/Gini/Gini_index_duplicates_{gene}.png', dpi=500)
    plt.close()
    print(f"file Gini_index_nodes_{gene}.png saved!")
def gini_run_edges_whole():
    print(f"Running with input: whole")
    gini_sjogren = []
    gini_healthy = []
    for file in os.listdir(sjogren_path):
        df = preprocessing(pd.read_csv(os.path.join(sjogren_path, file), sep='\t'))()
        cl = cluster(df)
        _, temp = cl.clonal_expansion_index_by_edges()
        gini_sjogren.append(temp)
    for file in os.listdir(healthy_path):
        df = preprocessing(pd.read_csv(os.path.join(healthy_path, file), sep='\t'))()
        cl = cluster(df)
        _, temp = cl.clonal_expansion_index_by_edges()
        gini_healthy.append(temp)

    max_length = max(len(gini_sjogren), len(gini_healthy))

    if len(gini_sjogren) < max_length:
        gini_sjogren.extend([np.nan] * (max_length - len(gini_sjogren)))
    if len(gini_healthy) < max_length:
        gini_healthy.extend([np.nan] * (max_length - len(gini_healthy)))

    df_gini = pd.DataFrame({'Healthy': gini_healthy, 'Sjogren': gini_sjogren})
    plt.figure(figsize=(4, 4))
    sns.set(style="whitegrid")
    ax = sns.boxplot(data=df_gini, palette=['lightgrey', 'skyblue'], width=0.35)
    ax.set(ylim=(0, 1), xlim=(0, 1))
    sns.stripplot(data=df_gini, jitter=True, color=".3")
    statannot.add_stat_annotation(ax, data=df_gini, box_pairs=[('Healthy', 'Sjogren')], test='t-test_ind',
                                  text_format='star', verbose=2, loc="inside")

    # Customize the plot
    plt.ylabel('Gini Index (edges)')
    plt.title('Whole')
    plt.tight_layout()
    plt.savefig(f'/home/scientist617/Grad_Proj/Data/Figures/Gini/Gini_index_nodes_whole.png', dpi=500)
    plt.close()
    print(f"file Gini_index_nodes_whole.png saved!")
def gini_run_duplicates_whole():
    print(f"Running with input: whole")
    gini_sjogren = []
    gini_healthy = []
    for file in os.listdir(sjogren_path):
        df = preprocessing(pd.read_csv(os.path.join(sjogren_path, file), sep='\t'))()
        cl = cluster(df)
        if len(df[df['v_call'] == gene]) == 0:
            continue
        _, temp = cl.clonal_expansion_index_by_duplicates()
        gini_sjogren.append(temp)
    for file in os.listdir(healthy_path):
        df = preprocessing(pd.read_csv(os.path.join(healthy_path, file), sep='\t'))()
        cl = cluster(df)
        if len(df[df['v_call'] == gene]) == 0:
            continue
        _, temp = cl.clonal_expansion_index_by_duplicates()
        gini_healthy.append(temp)

    max_length = max(len(gini_sjogren), len(gini_healthy))

    if len(gini_sjogren) < max_length:
        gini_sjogren.extend([np.nan] * (max_length - len(gini_sjogren)))
    if len(gini_healthy) < max_length:
        gini_healthy.extend([np.nan] * (max_length - len(gini_healthy)))

    df_gini = pd.DataFrame({'Healthy': gini_healthy, 'Sjogren': gini_sjogren})
    plt.figure(figsize=(4, 4))
    sns.set(style="whitegrid")
    ax = sns.boxplot(data=df_gini, palette=['lightgrey', 'skyblue'], width=0.35)
    ax.set(ylim=(0, 1), xlim=(0, 1))
    sns.stripplot(data=df_gini, jitter=True, color=".3")
    statannot.add_stat_annotation(ax, data=df_gini, box_pairs=[('Healthy', 'Sjogren')], test='t-test_ind',
                                  text_format='star', verbose=2, loc="inside")

    # Customize the plot
    plt.ylabel('Gini Index (duplicates)')
    plt.title(f'Whole')
    plt.tight_layout()
    plt.savefig(f'/home/scientist617/Grad_Proj/Data/Figures/Gini/Gini_index_duplicates_whole.png', dpi=500)
    plt.close()
    print(f"file Gini_index_nodes_whole.png saved!")
if __name__=="__main__":
    inputs = gene_list_all
    processes = []

    # for input_data in inputs:
    #     # process = Process(target=gini_run_edges, args=(input_data,))
    #     # processes.append(process)
    #     # process.start()
    #
    #     # process = Process(target=gini_run_duplicates, args=(input_data,))
    #     # processes.append(process)
    #     # process.start()

    process = Process(target=gini_run_edges_whole)
    process.start()
    processes.append(process)

    # process = Process(target=gini_run_duplicates_whole)
    # processes.append(process)
    # process.start()

    for process in processes:
        process.join()