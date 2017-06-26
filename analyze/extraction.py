import csv
import glob
import math
from pathlib import Path
import argparse

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from graph_tool.all import *


class ParameterExtraction:
    def __init__(self, graphs_path: str, export_file: str) -> None:
        super().__init__()
        self.graphs_path = Path(graphs_path)
        self.export_file = export_file

    def _iterate_edge_lists(self):
        if self.graphs_path.is_dir():
            for filename in glob.iglob(str(self.graphs_path) + '/**/*.edges', recursive=True):
                yield filename

    def _get_num_edges(self, graph: Graph):
        return graph.num_edges()

    def _get_num_vertices(self, graph: Graph):
        return graph.num_vertices()

    def _get_global_clustering_coefficient(self, graph: Graph):
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
        degree_sequence = sorted(nx.degree(G).values(), reverse=True)

    """
    DEGREE
    """

    def _get_avg_degree(self, graph: Graph):
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

    def run(self, do_non_np=True):
        columns = [self._get_num_edges,
                   self._get_num_vertices,
                   self._get_global_clustering_coefficient,
                   self._get_avg_degree,
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
        ]

        columns_nx_non_np = [
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
                print('Did graph {}'.format(file))
                column_values = [c(g) for c in columns]
                column_values += [c(gx) for c in columns_nx]
                if do_non_np:
                    column_values += [c(gx) for c in columns_nx_non_np]
                f.writerow(column_values)
            else:
                print('Graph was not connected, skipping ... {}'.format(file))
        f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extraction of diversity of graph properties')
    parser.add_argument('data_input', metavar='I', type=str,
                        help='Input data directory')
    parser.add_argument('data_output', metavar='O', type=str,
                        help='Output CSV file')
    parser.add_argument('-np',
                        help='Do graph properties algorithm that are not in NP')

    args = parser.parse_args()

    p = ParameterExtraction(args.data_input, args.data_output)

    if args.np:
        print('NP true')
        p.run(do_non_np=True)
    else:
        print('NP false')
        p.run(do_non_np=False)
