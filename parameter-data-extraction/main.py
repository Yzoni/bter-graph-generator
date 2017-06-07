from pathlib import Path
import glob
import csv
import matplotlib.pyplot as plt
from graph_tool.all import *
from sklearn import linear_model
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

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

    def _get_avg_deg(self, graph: Graph):
        return 2 * (graph.num_edges() / graph.num_vertices())

    def _get_pseudo_diameter(self, graph: Graph):
        return pseudo_diameter(graph)[0]

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
                   self._get_pseudo_diameter]

        f = csv.writer(open(self.export_file, 'w'))
        f.writerow(map(lambda x: x.__name__, columns))
        for file in self._iterate_edge_lists():
            with open(file, 'r') as c:
                dialect = csv.Sniffer().sniff(c.read(1024))
                csv_options = {"delimter": dialect.delimiter, "quatechar": dialect.quotechar}

            g = load_graph_from_csv(file, directed=False, csv_options=csv_options)
            column_values = [c(g) for c in columns]
            f.writerow(column_values)


def learn(csv_file: str):
    reg = linear_model.LinearRegression()
    arr = np.genfromtxt(csv_file, delimiter=',', skip_header=True).T
    y = np.vstack((arr[2], arr[3]))
    x = np.vstack((arr[4])).reshape((1114,))

    reg.fit(x.reshape(-1, 1), y.T)

    predict_size = 12
    predicted = [reg.predict(x) for x in np.arange(predict_size)]
    predicted = np.concatenate(predicted).T

    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(y[0], y[1], x)
    ax.plot(predicted[0], predicted[1], np.arange(predict_size), color='r')
    plt.show()


if __name__ == '__main__':
    p = ParameterExtraction('data', 'parameters.csv')
    p.run()
    learn('parameters.csv')
