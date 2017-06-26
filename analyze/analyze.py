import csv

import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate
from scipy.stats import spearmanr
from sklearn import linear_model

data = ['_get_num_edges', '_get_num_vertices', '_get_global_clustering_coefficient', '_get_avg_degree', '_get_pseudo_diameter', '_get_assortativity',
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

    print('Global clustering coefficient')
    cc = [[parameters[i]] + list(spearmanr(arr[2], arr[i])) for i in range(3, arr.shape[0])]
    print(tabulate(cc, headers=['parameter', 'rho', 'pval']))

    print('\n')

    print('Average degree distribution')
    avg_deg = [[parameters[i]] + list(spearmanr(arr[3], arr[i])) for i in range(4, arr.shape[0])]
    print(tabulate(avg_deg, headers=['parameter', 'rho', 'pval']))


def learn_diameter(csv_file: str):
    reg = linear_model.LinearRegression()
    arr = np.genfromtxt(csv_file, delimiter=',', skip_header=True).T
    y = np.vstack((arr[2], arr[3]))
    x = np.vstack((arr[4])).reshape((len(arr[4]),))

    reg.fit(x.reshape(-1, 1), y.T)

    predict_size = 20
    predicted = [reg.predict(x) for x in np.arange(predict_size)]
    predicted = np.concatenate(predicted).T

    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(y[0], y[1], x)
    ax.plot(predicted[0], predicted[1], np.arange(predict_size), color='r')

    ax.set_xlabel('Global clustering coefficient')
    ax.set_ylabel('Average global degree')
    ax.set_zlabel('Diameter')
    plt.show()

    return reg.coef_


def learn_scalar_assortativity(csv_file: str):
    reg = linear_model.LinearRegression()
    arr = np.genfromtxt(csv_file, delimiter=',', skip_header=True).T
    y = np.vstack((arr[2], arr[3]))
    x = np.vstack((arr[6])).reshape((len(arr[6]),))

    reg.fit(x.reshape(-1, 1), y.T)

    predicted = [reg.predict(x) for x in np.arange(-0.5, 0.5, 0.001)]
    predicted = np.concatenate(predicted).T

    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(y[0], y[1], x)
    ax.plot(predicted[0], predicted[1], np.arange(-0.5, 0.5, 0.001), color='r')

    ax.set_xlabel('Clustering coefficient')
    ax.set_ylabel('Average global degree')
    ax.set_zlabel('Scalar Assortativity')
    plt.show()

    return reg.coef_


def poly_fit_density(csv_file: str):
    arr = np.genfromtxt(csv_file, delimiter=',', names=True).T

    y = arr['_get_global_clustering_coefficient']
    x = arr['_get_density']
    xp = np.poly1d(np.polyfit(x, y, 2))
    yp = np.linspace(0, 1, 20)

    plt.scatter(y, x)
    plt.plot(xp(yp), yp, '-', color='r')
    plt.show()

def poly_fit_shortest_path(csv_file: str):
    arr = np.genfromtxt(csv_file, delimiter=',', names=True).T

    y = arr['_get_global_clustering_coefficient']
    x = arr['_get_average_shortest_path']
    xp = np.poly1d(np.polyfit(x, y, 2))
    yp = np.linspace(1, 2.6, 20)

    plt.scatter(y, x)
    plt.plot(xp(yp), yp, '-', color='r')
    plt.show()


def get_outliers(csv_file: str, column: str, cutoff: int):
    arr = np.genfromtxt(csv_file, delimiter=',', names=True).T
    return [{'index': index, 'value': x} for index, x in enumerate(arr[column]) if x < cutoff]


def plot_correlation(csv_file: str):
    arr = np.genfromtxt(csv_file, delimiter=',', names=True).T

    count = 0
    for column in data:
        if column is '_get_global_clustering_coefficient':
            continue
        count += 1
        x = arr['_get_global_clustering_coefficient']
        y = arr[column]

        plt.figure(1)
        plt.subplot(6, 3, count)
        plt.ylabel(column[5:])
        plt.scatter(x, y)
    plt.tight_layout()
    plt.show()

    count = 0
    for column in data:
        if column is '_get_avg_degree':
            continue
        count += 1

        x = arr['_get_avg_degree']
        y = arr[column]

        plt.figure(1)
        plt.subplot(6, 3, count)
        plt.ylabel(column[5:])
        plt.scatter(x, y)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    # correlation_spearman('parameters_proper.csv')
    plot_correlation('parameters_proper.csv')
    # poly_fit_density('parameters_proper.csv')
    # poly_fit_shortest_path('parameters_proper.csv')
    # print(get_outliers('parameters_proper.csv', '_get_ccd', 0.1))
