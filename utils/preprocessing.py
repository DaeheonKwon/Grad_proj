import pandas as pd
from collections import defaultdict


class preprocessing:

    def __init__(self, dataframe):
        self.dataframe = dataframe.copy()
        self.trim_alleles()
        self.grouped_dataframe = self.dataframe.groupby(['v_call', 'j_call', 'cdr3_aa'])

    def __call__(self):
        return self.dataframe

    def trim_alleles(self) -> object:
        self.dataframe['v_call'] = self.dataframe['v_call'].apply(lambda x: x.split('*')[0])
        self.dataframe['j_call'] = self.dataframe['j_call'].apply(lambda x: x.split('*')[0])
        self.dataframe['c_call'] = self.dataframe['c_call'].apply(lambda x: x.split('*')[0])

    def get_naive(self):
        naive = self.dataframe[self.dataframe['v_alignment_mutation'] < 2]
        naive = naive[naive['c_call'].isin(['IGHM', 'IGHD'])]
        return naive

    def get_mature(self):
        mature = self.dataframe[self.dataframe['v_alignment_mutation'] >= 2]
        # consider cases that c_call is other than M/D, but SHM is under 2


        return mature