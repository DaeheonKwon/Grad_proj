import pandas as pd
from collections import defaultdict


class cluster:

    def __init__(self, dataframe):
        self.dataframe = dataframe
        self.dataframe['cdr3_aa_len'] = self.dataframe['cdr3_aa'].apply(lambda x: len(x))
        self.v_call = self.dataframe['v_call']
        self.j_call = self.dataframe['j_call']
        self.c_call = self.dataframe['c_call']
        self.trim_alleles()
        self.grouped_dataframe = self.dataframe.groupby(['v_call', 'j_call', 'cdr3_aa_len'])
        self.grouped_dataframe = list(self.grouped_dataframe)

    def __call__(self):
        grouped_hamming_distance_map = {}
        for i in range(len(self.grouped_dataframe)):
            grouped_hamming_distance_map[i] = self.search_adjacency(self.grouped_dataframe[i])
        return grouped_hamming_distance_map    

    def trim_alleles(self) -> object:
        for i in range(len(self.v_call)):
            for j in range(len(self.v_call[i])):
                if self.v_call[i][j] == '*':
                    self.v_call[i] = self.v_call[i][:j]
                    break
        for i in range(len(self.j_call)):
            for j in range(len(self.j_call[i])):
                if self.j_call[i][j] == '*':
                    self.j_call[i] = self.j_call[i][:j]
                    break
        for i in range(len(self.c_call)):
            for j in range(len(self.c_call[i])):
                if self.c_call[i][j] == '*':
                    self.c_call[i] = self.c_call[i][:j]
                    break

    def search_adjacency(self, grouped_dataframe_subgroup):
        sequence = grouped_dataframe_subgroup[1]['cdr3_aa']
        sequence_length: object = grouped_dataframe_subgroup[1]['cdr3_aa_len'].iloc[0]
        id = grouped_dataframe_subgroup[1]['sequence_id']
        hamming_distance_map = defaultdict(list)

        for i in range(len(sequence)):
            for j in range(i + 1, len(sequence)):
                hamming_distance = 0
                for k in range(sequence_length):
                    if sequence.iloc[i][k] != sequence.iloc[j][k]:
                        hamming_distance += 1
                    if hamming_distance > 1:
                        break
                if hamming_distance == 1:
                    hamming_distance_map[id.iloc[i]].append(id.iloc[j])
                    hamming_distance_map[id.iloc[j]].append(id.iloc[i])
        return hamming_distance_map
