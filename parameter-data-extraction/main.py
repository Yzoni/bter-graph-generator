from pathlib import Path
import glob
import csv
import matplotlib.pyplot as plt
from graph_tool.all import *
from sklearn import linear_model
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import math
import networkx as nx
from scipy.stats import spearmanr
from tabulate import tabulate

class ParameterExtraction:
    def __init__(self, graphs_path: str, export_file: str) -> None:
        super().__init__()
        self.graphs_path = Path(graphs_path)
        self.export_file = export_file

    def _iterate_edge_lists(self):
        if self.graphs_path.is_dir():
            for filename in glob.iglob(str(self.graphs_path) + '/**/*.edges', recursive=True):
                print(filename)
                yield filename

    def _get_num_edges(self, graph: Graph):
        return graph.num_edges()

    def _get_num_vertices(self, graph: Graph):
        return graph.num_vertices()

    def _get_ccd(self, graph: Graph):
        return global_clustering(graph)[0]

    def _get_pseudo_diameter(self, graph: Graph):
        return pseudo_diameter(graph)[0]

    def _get_assortativity(self, graph: Graph):
        return assortativity(graph, 'total')[0]

    def _get_scalar_assortativity(self, graph: Graph):
        return scalar_assortativity(graph, 'total')[0]

    def _get_density(self, graph: Graph):
        return (2 * graph.num_edges()) / (graph.num_vertices() * (graph.num_vertices() - 1))

    def _get_average_shortest_path(self, graph: nx.Graph):
        return nx.average_shortest_path_length(graph)

    # Degree distribution log normal?
    def get_degree_sequence(self, graph: nx.Graph):
        degree_sequence=sorted(nx.degree(G).values(),reverse=True)

    """
    DEGREE
    """
    def _get_avg_deg(self, graph: Graph):
        return 2 * (graph.num_edges() / graph.num_vertices())

    def _get_min_deg(self, graph: Graph):
        min_degree = math.inf
        for v in graph.vertices():
            degree = v.out_degree()
            if (degree < min_degree):
                min_degree = degree
        return min_degree

    def _get_max_deg(self, graph: Graph):
        max_degree = -1
        for v in graph.vertices():
            degree = v.out_degree()
            if (degree > max_degree):
                max_degree = degree
        return max_degree

    """
    VERTEX BETWEENNESS
    """
    def _get_max_vertex_betweenness(self, graph: Graph):
        vertex_betweenness, edge_betweenness = betweenness(graph)
        return max(vertex_betweenness.get_array())

    def _get_min_vertex_betweenness(self, graph: Graph):
        vertex_betweenness, edge_betweenness = betweenness(graph)
        return min(vertex_betweenness.get_array())

    def _get_mean_vertex_betweenness(self, graph: Graph):
        vertex_betweenness, edge_betweenness = betweenness(graph)
        return np.mean(vertex_betweenness.get_array())

    def _get_std_vertex_betweenness(self, graph: Graph):
        vertex_betweenness, edge_betweenness = betweenness(graph)
        return np.std(vertex_betweenness.get_array())

    """
    EDGE BETWEENNESS
    """
    def _get_max_edge_betweenness(self, graph: Graph):
        vertex_betweenness, edge_betweenness = betweenness(graph)
        return max(edge_betweenness.get_array())

    def _get_min_edge_betweenness(self, graph: Graph):
        vertex_betweenness, edge_betweenness = betweenness(graph)
        return min(edge_betweenness.get_array())

    def _get_mean_edge_betweenness(self, graph: Graph):
        vertex_betweenness, edge_betweenness = betweenness(graph)
        return np.mean(edge_betweenness.get_array())

    def _get_std_edge_betweenness(self, graph: Graph):
        vertex_betweenness, edge_betweenness = betweenness(graph)
        return np.std(edge_betweenness.get_array())

    """
    CLIQUE
    """
    def _get_clique_number(self, graph: nx.Graph):
        return nx.graph_clique_number(graph)

    def _get_number_of_cliques(self, graph: nx.Graph):
        return nx.graph_number_of_cliques(graph)


    """
    PLOT
    """
    def plot_data(self, csv_file: str):
        p = Path(csv_file)
        if p.exists():
            with open(str(p)) as csv_file:
                reader = csv.reader(csv_file, delimiter=',')
                next(reader)
                ccd = list()
                avd = list()
                for row in reader:
                    ccd.append(row[2])
                    avd.append(row[3])
            plt.scatter(ccd, avd)
            plt.show()

    def run(self):
        columns = [self._get_num_edges,
                   self._get_num_vertices,
                   self._get_ccd,
                   self._get_avg_deg,
                   self._get_pseudo_diameter,
                   self._get_assortativity,
                   self._get_scalar_assortativity,
                   self._get_density,
                   self._get_min_deg,
                   self._get_max_deg,
                   self._get_max_vertex_betweenness,
                   self._get_min_vertex_betweenness,
                   self._get_mean_vertex_betweenness,
                   self._get_std_vertex_betweenness,
                   self._get_max_edge_betweenness,
                   self._get_min_edge_betweenness,
                   self._get_mean_edge_betweenness,
                   self._get_std_edge_betweenness]

        columns_nx = [
                    self._get_clique_number,
                    self._get_number_of_cliques,
                    self._get_average_shortest_path
        ]

        f = csv.writer(open(self.export_file, 'w'))
        f.writerow(map(lambda x: x.__name__, columns + columns_nx))
        for file in self._iterate_edge_lists():

            # Get data properties
            with open(file, 'r') as c:
                dialect = csv.Sniffer().sniff(c.read(1024))
                csv_options = {'delimiter': dialect.delimiter, 'quotechar': dialect.quotechar}

            # Create graph and write specified columns to new csv file
            g = load_graph_from_csv(file, directed=False, csv_options=csv_options)
            gx = nx.read_edgelist(file, delimiter=csv_options['delimiter'], create_using=nx.Graph(), nodetype=int)
            if nx.is_connected(gx):
                column_values = [c(g) for c in columns]
                column_values += [c(gx) for c in columns_nx]
                f.writerow(column_values)
        f.close()

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

    ax.set_xlabel('Clustering coefficient')
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

    y = arr['_get_ccd']
    x = arr['_get_density']

    plt.scatter(x, y)
    plt.show()

if __name__ == '__main__':
    p = ParameterExtraction('data', 'parameters.csv')
    # p.run()
    # learn_diameter('parameters.csv')
    # learn_scalar_assortativity('parameters.csv')
    correlation_spearman('parameters.csv')
