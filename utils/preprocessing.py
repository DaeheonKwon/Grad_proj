import pandas as pd
import numpy as np
from collections import defaultdict


class preprocessing:

    def __init__(self, dataframe):
        self.dataframe = dataframe.copy()
        self.dataframe.fillna('*', inplace=True)
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
        mature_1 = self.dataframe[self.dataframe['v_alignment_mutation'] >= 2]
        mature_2 = self.dataframe[self.dataframe['v_alignment_mutation'] < 2]
        mature_2 = mature_2[~mature_2['c_call'].isin(['IGHM', 'IGHD'])]
        return pd.concat([mature_1, mature_2])

    def get_Vgene(self, Vgene):
        return self.dataframe[self.dataframe['v_call'] == Vgene]

    def get_mature_usage(self):
        v_counts = self.get_mature()['v_call'].value_counts()
        j_counts = self.get_mature()['j_call'].value_counts()
        c_counts = self.get_mature()['c_call'].value_counts()

        return v_counts, j_counts, c_counts

    def get_naive_usage(self):
        v_counts = self.get_naive()['v_call'].value_counts()
        j_counts = self.get_naive()['j_call'].value_counts()
        c_counts = self.get_naive()['c_call'].value_counts()

        return v_counts, j_counts, c_counts
