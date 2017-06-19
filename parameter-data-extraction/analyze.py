import csv

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from tabulate import tabulate

data = ['_get_num_edges', '_get_num_vertices', '_get_ccd', '_get_avg_deg', '_get_pseudo_diameter', '_get_assortativity',
        '_get_scalar_assortativity', '_get_density', '_get_min_deg', '_get_max_deg', '_get_max_vertex_betweenness',
        '_get_min_vertex_betweenness', '_get_mean_vertex_betweenness', '_get_std_vertex_betweenness',
        '_get_max_edge_betweenness', '_get_min_edge_betweenness', '_get_mean_edge_betweenness',
        '_get_std_edge_betweenness', '_get_average_shortest_path']


def correlation_spearman(csv_file: str):
    """
    Generates a table of spearman correlations for each column of csv

    :param csv_file: csv file name path
    :return:
    """
    with open(csv_file) as f:
        reader = csv.reader(f)
        parameters = next(reader)

    arr = np.genfromtxt(csv_file, delimiter=',', skip_header=True).T

    print('Clustering coefficient')
    cc = [[parameters[i]] + list(spearmanr(arr[2], arr[i])) for i in range(3, arr.shape[0])]
    print(tabulate(cc, headers=['parameter', 'rho', 'pval']))

    print('\n')

    print('Average degree distribution')
    avg_deg = [[parameters[i]] + list(spearmanr(arr[3], arr[i])) for i in range(4, arr.shape[0])]
    print(tabulate(avg_deg, headers=['parameter', 'rho', 'pval']))


def plot_correlation(csv_file: str):
    arr = np.genfromtxt(csv_file, delimiter=',', names=True).T

    for idx, column in enumerate(data):
        y = arr['_get_ccd']
        x = arr[column]

        plt.figure(1)
        plt.subplot(5, 4, idx + 1)
        plt.title(column)
        plt.scatter(x, y)
    plt.show()

    for idx, column in enumerate(data):
        y = arr['_get_avg_deg']
        x = arr[column]

        plt.figure(1)
        plt.subplot(5, 4, idx + 1)
        plt.title(column)
        plt.scatter(x, y)
    plt.show()


if __name__ == '__main__':
    correlation_spearman('parameters_proper.csv')
    plot_correlation('parameters_proper.csv')
