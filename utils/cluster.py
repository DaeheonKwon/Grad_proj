import pandas as pd
from collections import defaultdict


class cluster:

    def __init__(self, dataframe):
        self.dataframe = dataframe
        self.dataframe['cdr3_aa_len'] = self.dataframe['cdr3_aa'].apply(lambda x: len(x))
        self.trim_alleles()
        self.grouped_dataframe = self.dataframe.groupby(['v_call', 'j_call', 'cdr3_aa_len'])
        self.grouped_dataframe = list(self.grouped_dataframe)

    def __call__(self):
        grouped_hamming_distance_map = {}
        for i in range(len(self.grouped_dataframe)):
            grouped_hamming_distance_map[i] = self.search_adjacency(self.grouped_dataframe[i])
        return grouped_hamming_distance_map    

    def trim_alleles(self) -> object:
        self.dataframe['v_call'] = self.dataframe['v_call'].apply(lambda x: x.split('*')[0])
        self.dataframe['j_call'] = self.dataframe['j_call'].apply(lambda x: x.split('*')[0])
        self.dataframe['c_call'] = self.dataframe['c_call'].apply(lambda x: x.split('*')[0])

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
