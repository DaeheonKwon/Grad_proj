import pandas as pd
from itertools import combinations


class recombinations:
    def __init__(self, dataframe):
        self.dataframe = dataframe
        self.v_call = self.dataframe['v_call']
        self.j_call = self.dataframe['j_call']
        self.c_call = self.dataframe['c_call']
        self.trim_alleles()

    def __call__(self):
        recombination_types = {}

        for _, constant_genes in self.check_for_CSR():
            unique_constant_genes = set(constant_genes)

            for combo in combinations(unique_constant_genes, 2):
                if 'IGHGP' in combo:
                    continue
                if ('IGHA1' in combo and 'IGHA2' in combo) or ('IGHM' in combo and 'IGHD' in combo) or (
                        'IGHG1' in combo and 'IGHG2' in combo):
                    continue

                recombination_type = "-".join(sorted(combo))
                recombination_types[recombination_type] = recombination_types.get(recombination_type, 0) + 1

        return recombination_types

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

    def check_for_CSR(self):
        clones = self.dataframe.groupby(['v_call', 'j_call', 'cdr3_aa'])
        clones = list(clones)
        clones_with_CSR = []

        for i in range(len(clones)):
            isCSR, Cgenes = self.check_multiple_Cgene(clones[i][1])
            if isCSR:
                clones_with_CSR.append((clones[i][0], Cgenes))
        return clones_with_CSR

    def check_multiple_Cgene(self, clone_dataframe):
        isCSR = clone_dataframe['c_call'].nunique() > 1
        Cgenes = clone_dataframe['c_call'].unique()

        return isCSR, Cgenes
