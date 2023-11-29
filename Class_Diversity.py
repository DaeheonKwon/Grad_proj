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
healthy_path = '/home/scientist617/Grad_Proj/Data/new_in_house_healthy/'


# %%
def checkCDRpos(respos, chain):
    for cdr, pos in cdrpos[chain].items():
        if compare_strings(pos[0], respos) & compare_strings(respos, pos[1]):
            return True, cdr
    return False, None


def inVgene(respos):
    if compare_strings(respos, '113'):
        return True
    return False


def compare_strings(str1, str2):
    # Extract numerical values from the input strings using regular expressions
    num1 = int(re.search(r'\d+', str1).group())
    num2 = int(re.search(r'\d+', str2).group())

    # If the numerical values are equal, compare the remaining parts of the strings
    if num1 == num2:
        alpha1 = re.search(r'[A-Za-z]+', str1).group() if re.search(r'[A-Za-z]+', str1) else ""
        alpha2 = re.search(r'[A-Za-z]+', str2).group() if re.search(r'[A-Za-z]+', str2) else ""

        # Compare the alphabetic parts using string comparison
        return alpha1 <= alpha2

    # Compare the numerical values
    return num1 <= num2


def custom_sort(strings):
    def extract_key(s):
        # Split the string into numeric and non-numeric parts
        parts = re.split(r'(\d+)', s)

        # Convert numeric parts to integers, leave non-numeric parts as they are
        converted_parts = [int(part) if part.isdigit() else part for part in parts]

        return converted_parts

    return sorted(strings, key=extract_key)


def concat_usage(path):
    data = defaultdict(pd.DataFrame)
    for file in os.listdir(path):
        df = preprocessing(pd.read_csv(os.path.join(path, file), sep='\t'))
        v, _, _ = df.get_mature_usage()
        v1, _, _ = df.get_naive_usage()
        mature_data = v.div(v.sum())
        naive_data = v1.div(v1.sum())
        # Combine the data into a single DataFrame
        df = pd.DataFrame({'Naive Repertoire': naive_data, 'Mature Repertoire': mature_data})
        data[file.split('_')[0]] = df
    dataframes = []
    for patient, repertoires in data.items():
        dataframes.append(repertoires)
    combined_df = pd.concat(dataframes)
    return combined_df


def convert_pvalue_to_asterisks(pvalue):
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
    return "ns"


def is_pattern_present(input_string, pattern):
    # Replace "_" in the pattern with "." in the regex to match any character
    pattern = str(pattern).replace("_", ".")
    # Create a regex pattern
    regex_pattern = re.compile(pattern)
    # Search for the pattern in the input string
    match = regex_pattern.search(str(input_string))

    return match is not None


def levenshtein_distance(str1, str2):
    str1 = str(str1)
    str2 = str(str2)
    len_str1 = len(str1) + 1
    len_str2 = len(str2) + 1

    # Initialize matrix to store distances
    distance_matrix = [[0 for _ in range(len_str2)] for _ in range(len_str1)]

    # Initialize the first row and column
    for i in range(len_str1):
        distance_matrix[i][0] = i
    for j in range(len_str2):
        distance_matrix[0][j] = j

    # Fill in the matrix
    for i in range(1, len_str1):
        for j in range(1, len_str2):
            cost = 0 if str1[i - 1] == str2[j - 1] else 1

            distance_matrix[i][j] = min(
                distance_matrix[i - 1][j] + 1,  # Deletion
                distance_matrix[i][j - 1] + 1,  # Insertion
                distance_matrix[i - 1][j - 1] + cost  # Substitution
            )

    return distance_matrix[len_str1 - 1][len_str2 - 1] - abs(len_str1 - len_str2), abs(len_str1 - len_str2)


def sequence_similarity(seq1, seq2):
    sub, indel = levenshtein_distance(seq1, seq2)
    sequence_similarity = pow(0.3, sub) * pow(0.1, indel)
    return sequence_similarity


# %%
path = '/home/scientist617/Grad_Proj/Data/PDB structures/vh_match_100/'
df = pd.read_csv(os.path.join(path, 'vh_match_filtered_immuneSys.tsv'), sep='\t')
df = df.groupby('pdb').first()
pdb_parser = PDBParser(QUIET=True)
structure = pdb_parser.get_structure('5xmh', os.path.join(path, 'chothia', '5xmh.pdb'))
query_string_whole = ""
query_string_CDR1 = ''
query_string_FR2 = ''
query_string_CDR2 = ''
query_string_FR3 = ''
query_string_CDR3 = ''
for i in structure[0][df.loc['5xmh', 'Hchain']]:
    if 25 < i.get_id()[1] & i.get_id()[1] < 102:
        query_string_whole = query_string_whole + seq1(i.get_resname())
    if 25 < i.get_id()[1] & i.get_id()[1] < 33:
        query_string_CDR1 = query_string_CDR1 + seq1(i.get_resname())
    if 32 < i.get_id()[1] & i.get_id()[1] < 52:
        query_string_FR2 = query_string_FR2 + seq1(i.get_resname())
    if 51 < i.get_id()[1] & i.get_id()[1] < 57:
        query_string_CDR2 = query_string_CDR2 + seq1(i.get_resname())
    if 56 < i.get_id()[1] & i.get_id()[1] < 96:
        query_string_FR3 = query_string_FR3 + seq1(i.get_resname())
    if 95 < i.get_id()[1] & i.get_id()[1] < 102:
        query_string_CDR3 = query_string_CDR3 + seq1(i.get_resname())
# %%
motifs = [
    # Seq('GTF'),
    # Seq('GGI'),
    # Seq('PGQ'),
    # Seq('GII'),
    # Seq('AGT'),
    # Seq('SAG'),
    # Seq('GTP'),
    # Seq('W_G'),
    # Seq('Q__E'),
    # Seq(query_string_CDR1),
    # Seq(query_string_FR2),
    # Seq(query_string_CDR2),
    # Seq(query_string_FR3),
    Seq(query_string_CDR3),
    # Seq(query_string_whole)
]

motif_types = ['CDR3']


def get_motif_score(path, motif, type, q):
    output_score = []

    for file in os.listdir(path):
        if os.path.isdir(os.path.join(path, file)):
            continue
        pre = preprocessing(pd.read_csv(os.path.join(path, file), sep='\t'))
        df_naive = pre.get_naive()
        df_mature = pre.get_mature()
        # trim CDRs to match chothia numbering system

        def update_score(df):
            similarity_score = 0
            similarity_score_q = 1

            for i, row in df.iterrows():
                frequency = row['duplicate_count']/df['duplicate_count'].sum(axis=0)

                str_pos = []
                seq = Seq(row['sequence_aa'])
                seq_1 = Seq(row['cdr1_aa'][:-1])
                seq_2 = Seq(row['cdr2_aa'][1:-1])
                seq_3 = Seq(row['cdr3_aa'][:-1])
                cdr_pos = [seq.find(seq_1), seq.find(seq_2), seq.find(seq_3)]
                cdr_length = [len(seq_1), len(seq_2), len(seq[seq.find(seq_2) + len(seq_2): seq.find(seq_3)]), len(seq_3)]
                cdr_length_ori = [int(cdrpos['H']['H1'][1]) - int(cdrpos['H']['H1'][0]) + 1,
                                  int(cdrpos['H']['H2'][1]) - int(cdrpos['H']['H2'][0]) + 1,
                                  int(cdrpos['H']['H3'][0]) - int(cdrpos['H']['H2'][1]) - 1,
                                  int(cdrpos['H']['H3'][1]) - int(cdrpos['H']['H3'][0]) + 1]

                chothia_num = list(range(26 - seq.find(seq_1), 102 + len(seq[seq.find(seq_3) + len(seq_3):])))
                chothia_num = [str(i) for i in chothia_num]
                duplicate_chothia = ['31', '52', '82', '100']
                for i in range(len(cdr_length)):
                    for j in range(cdr_length[i] - cdr_length_ori[i]):
                        chothia_num.append(duplicate_chothia[i] + chr(ord('A') + j))
                chothia_num = custom_sort(chothia_num)
                chothia_num_seq = list(zip(seq, chothia_num))
                for i, (aa, num) in enumerate(chothia_num_seq):
                    if num == '26':
                        str_pos.append(i)
                    if num == '32':
                        str_pos.append(i)
                    if num == '52':
                        str_pos.append(i)
                    if num == '56':
                        str_pos.append(i)
                    if num == '96':
                        str_pos.append(i)
                    if num == '101':
                        str_pos.append(i)

                if len(str_pos) < 6:
                    continue

                chothiaseq = ''.join([i[0] for i in chothia_num_seq])
                # search for XXX formatted motifs, and append sequence_id to list
                if type == 'CDR1':
                    if q == 1:
                        similarity_score_q = similarity_score_q * pow(sequence_similarity(
                            chothiaseq[str_pos[0]:str_pos[1] + 1], motif), frequency)
                    else:
                        similarity_score = similarity_score + frequency * pow(
                            sequence_similarity(chothiaseq[str_pos[0]:str_pos[1] + 1], motif), q - 1)
                elif type == 'FR2':
                    if q == 1:
                        similarity_score_q = similarity_score_q * pow(sequence_similarity(
                            chothiaseq[str_pos[1] + 1:str_pos[2]], motif), frequency)
                    else:
                        similarity_score = similarity_score + frequency * pow(
                            sequence_similarity(chothiaseq[str_pos[1] + 1:str_pos[2]], motif), q - 1)
                elif type == 'CDR2':
                    if q == 1:
                        similarity_score_q = similarity_score_q * pow(sequence_similarity(
                            chothiaseq[str_pos[2]:str_pos[3] + 1], motif), frequency)
                    else:
                        similarity_score = similarity_score + frequency * pow(
                            sequence_similarity(chothiaseq[str_pos[2]:str_pos[3] + 1], motif), q - 1)
                elif type == 'FR3':
                    if q == 1:
                        similarity_score_q = similarity_score_q * pow(sequence_similarity(
                            chothiaseq[str_pos[3] + 1:str_pos[4]], motif), frequency)
                    else:
                        similarity_score = similarity_score + frequency * pow(
                            sequence_similarity(chothiaseq[str_pos[3] + 1:str_pos[4]], motif), q - 1)
                elif type == 'CDR3':
                    if q == 1:
                        similarity_score_q = similarity_score_q * pow(sequence_similarity(
                            chothiaseq[str_pos[4]:str_pos[5] + 1], motif), frequency)
                    else:
                        similarity_score = similarity_score + frequency * pow(
                            sequence_similarity(chothiaseq[str_pos[4]:str_pos[5] + 1], motif), q - 1)
                elif type == 'whole':
                    if q == 1:
                        similarity_score_q = similarity_score_q * pow(sequence_similarity(
                            chothiaseq[str_pos[0]:str_pos[5] + 1], motif), frequency)
                    else:
                        similarity_score = similarity_score + frequency * pow(
                            sequence_similarity(chothiaseq[str_pos[0]:str_pos[5] + 1], motif), q - 1)
            if q==1:
                return similarity_score_q
            else:
                return similarity_score

        if q == 1:
            output_score.append([1 / update_score(df_naive), 1 / update_score(df_mature),
                                 1 / (update_score(pre()))])
        else:
            output_score.append([pow(update_score(df_naive), 1 / (1 - q)), pow(update_score(df_mature), 1 / (1 - q)),
                                 pow(update_score(pre()), 1 / (1 - q))])
        print(f"{file}, {type} finished!")
    return output_score

def run_script(q):
    print(f"Running with input: {q}")
    pd.DataFrame([get_motif_score(sjogren_path, motif, motif_types[i], q) for i, motif in enumerate(motifs)]).to_csv(
        f'/home/scientist617/Grad_Proj/Data/Simm_Index_diff_q/sjogren_score_q_{q}.csv', index=False)
def run_script_ctrl(q):
    print(f'Running with input: {q}')
    pd.DataFrame([get_motif_score(healthy_path, motif, motif_types[i], q) for i, motif in enumerate(motifs)]).to_csv(
        f'/home/scientist617/Grad_Proj/Data/Simm_Index_diff_q/healthy_score_q_{q}.csv', index=False)


if __name__=="__main__":
    inputs = [20, 30, 40, 50, 60, 70, 80, 90, 100, 1000, 10000]
    processes = []

    for input_data in inputs:
        process = Process(target=run_script, args=(input_data,))
        processes.append(process)
        process.start()

    for input_data in inputs:
        process = Process(target=run_script_ctrl, args=(input_data,))
        processes.append(process)
        process.start()

    for process in processes:
        process.join()