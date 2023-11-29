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

sjogren_path = '/home/scientist617/Grad_Proj/Data/new_sjogren_file'
healthy_path = '/home/scientist617/Grad_Proj/Data/in_house_healthy'

def run_script(input_data):
    print(f"running with input: {input_data}")
    df_paths = []
    paths = [sjogren_path, healthy_path]
    for path in paths:
        files = os.listdir(path)
        series_list = []

        for file in files:
            if file.split('_')[0] == 'S8':
                continue
            df = preprocessing(pd.read_csv(os.path.join(path, file), sep='\t'))
            if input_data ==  'mature':
                df = df.get_mature()
            elif input_data == 'naive':
                df = df.get_naive()
            elif input_data == 'whole':
                df = df()
            series = pd.DataFrame(df.groupby('v_call')['cdr3_aa'].apply(lambda x: x.apply(lambda x: len(x)))).reset_index()

            for vgene in ['IGHV1-18', 'IGHV1-24','IGHV1-69', 'IGHV3-23', 'IGHV3-43','IGHV3-53', 'IGHV3-72', 'IGHV4-34', 'IGHV5-51', 'IGHV6-1']:
                for i in range(0, 100):
                    series_list.append(series[series['v_call'] == vgene].sample(20))

        temp = pd.DataFrame(pd.concat(series_list, axis=0)).reset_index()
        temp = temp[['v_call', 'cdr3_aa']]

        df_paths.append(temp)

    df_paths[0].reset_index(drop=False, inplace=True)
    df_paths[0]['index'] = ['Sjogren']*len(df_paths[0])
    df_paths[1].reset_index(drop=False, inplace=True)
    df_paths[1]['index'] = ['Healthy']*len(df_paths[1])
    df_paths[0] = df_paths[0][['index', 'v_call', 'cdr3_aa']]
    df_paths[1] = df_paths[1][['index', 'v_call', 'cdr3_aa']]
    merged_df = pd.concat(df_paths, axis=0)
    merged_df.columns = ['Group', 'V Genes', 'CDR3 Length']

    plt.figure(figsize=(10, 4))
    data = merged_df

    # Create the violin plot for both 'Healthy' and 'Sjogren' columns
    ax = sns.violinplot(data=data, x='V Genes', y='CDR3 Length', hue='Group', palette='Set3')

    box_pairs = []
    variables = data['V Genes'].unique()
    groups = data['Group'].unique()

    for var in variables:
        box_pairs.append(((var, groups[0]), (var, groups[1])))

    statannot.add_stat_annotation(ax, data=data, x='V Genes', y='CDR3 Length', hue='Group', box_pairs=box_pairs,
                                  test='Mann-Whitney', text_format='star', verbose=2, loc="outside")

    # plt.setp(ax.xaxis.get_majorticklabels(), rotation=-45, ha="left" )
    sns.move_legend(
        ax, "upper right", ncol=1, title=None, frameon=True
    )
    plt.tight_layout()
    sns.despine()
    plt.savefig(f'/home/scientist617/Grad_Proj/Data/Figures/CDR3_length_pooled/CDR3_length_{input_data}.png', dpi=500)

if __name__=="__main__":
    inputs = ['mature', 'naive', 'whole']
    processes = []

    for input_data in inputs:
        process = Process(target=run_script, args=(input_data,))
        processes.append(process)
        process.start()

    for process in processes:
        process.join()

