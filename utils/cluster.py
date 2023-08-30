import pandas as pd
from collections import defaultdict
import networkx as nx
import community
import numpy as np
import matplotlib.pyplot as plt

class cluster:

    def __init__(self, dataframe):
        self.dataframe = dataframe.copy()
        self.dataframe['cdr3_aa_len'] = self.dataframe['cdr3_aa'].apply(lambda x: len(x))
        self.trim_alleles()
        self.grouped_dataframe = self.dataframe.groupby(['v_call', 'j_call', 'cdr3_aa_len'])
        self.grouped_dataframe = list(self.grouped_dataframe)

    def trim_alleles(self) -> object:
        self.dataframe['v_call'] = self.dataframe['v_call'].apply(lambda x: x.split('*')[0])
        self.dataframe['j_call'] = self.dataframe['j_call'].apply(lambda x: x.split('*')[0])
        self.dataframe['c_call'] = self.dataframe['c_call'].apply(lambda x: x.split('*')[0])

    def distance_map(self):
        grouped_hamming_distance_map = {}
        for i in range(len(self.grouped_dataframe)):
            grouped_hamming_distance_map[i] = self.search_adjacency(self.grouped_dataframe[i])
        return grouped_hamming_distance_map

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

    def create_graph(self):
        G = nx.Graph()
        distance_map = self.distance_map()
        for i in range(len(distance_map)):
            for seq, neighbors in distance_map[i].items():
                G.add_node(seq)
                for neighbor in neighbors:
                    G.add_edge(seq, neighbor)
        return G

    def draw_graph(self):
        G = self.create_graph()
        layout = nx.spring_layout(G)
        nx.draw(G, layout, with_labels = False, node_size = 1, node_color = 'skyblue', edge_color = 'black', font_size = 1)
        plt.show()

    def clonal_expansion_index_by_edges(self):
        G = self.create_graph()
        edge_counts = [G.degree(vertex) for vertex in G.nodes]
        edge_counts.sort()

        cumulative_counts = np.cumsum(edge_counts)
        lorenz_curve_values = cumulative_counts / cumulative_counts[-1]
        area_lorenz_curve = np.trapz(lorenz_curve_values)
        area_perfect_equality = 0.5 * len(edge_counts)
        gini_coefficient = (area_perfect_equality - area_lorenz_curve) / area_perfect_equality

        return edge_counts, gini_coefficient

    def clonal_expansion_index_by_duplicates(self):
        sequence_counts = self.dataframe.groupby('sequence')['duplicate_count'].sum().tolist()
        sequence_counts.sort()

        cumulative_counts = np.cumsum(sequence_counts)
        lorenz_curve_values = cumulative_counts / cumulative_counts[-1]
        area_lorenz_curve = np.trapz(lorenz_curve_values)
        area_perfect_equality = 0.5 * len(sequence_counts)
        gini_coefficient = (area_perfect_equality - area_lorenz_curve) / area_perfect_equality

        return sequence_counts, gini_coefficient

    def clonal_diversification_index(self):
        G = self.create_graph()
        partition = community.best_partition(G)
        cluster_dict = {}
        for node, cluster_id in partition.items():
            if cluster_id not in cluster_dict:
                cluster_dict[cluster_id] = []
            cluster_dict[cluster_id].append(node)

        clusters = list(cluster_dict.values())

        cluster_size = []
        for cluster in clusters:
            cluster_size.append(len(cluster))
        cluster_size.sort()

        total_size = sum(cluster_size)
        proportions = [size/total_size for size in cluster_size]
        shannon_entropy = -sum([p * np.log2(p) for p in proportions])

        return cluster_size, shannon_entropy