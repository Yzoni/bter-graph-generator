import subprocess
import numpy as np
from pathlib import Path


def generate_analyze_graphs():
    number_of_vertices = 1000
    max_cluster_coefficient = 0.95

    range_max_degree_bound = np.arange(100, 250, 25)
    range_average_degree = np.arange(10, 100, 10)
    range_clustering_coefficient = np.arange(0.1, 0.9, 0.1)

    for db in range_max_degree_bound:
        for ad in range_average_degree:
            for ccf in range_clustering_coefficient:
                graph_out_dir = Path('generated_graphs')
                if not graph_out_dir.exists():
                    graph_out_dir.mkdir()

                graph_out_file = Path(str(db) + '_' + str(ad) + '_' + str(ccf) + '_graph.edges')
                args = ('../cmake-build-debug/bter_run',
                        '-n', str(number_of_vertices),
                        '-b', str(db),
                        '-d', str(ad),
                        '-m', str(max_cluster_coefficient),
                        '-c', str(ccf),
                        '--gpu',
                        '-f', str(graph_out_dir / graph_out_file))

                print(str(db) + '\t' + str(ad) + '\t', str(ccf))
                popen = subprocess.Popen(args, stdout=subprocess.PIPE)
                popen.wait()
                output = popen.stdout.read()
                print(output)


def test_graphs_from_real_graphs():
    number_of_vertices = 1000
    max_cluster_coefficient = 0.90
    max_degree_bound = 100
    average_degree = 5

    range_density = np.arange(0, 1.1, 0.1)
    range_shortest_path = np.arange(1, 2.5, 0.1)

    for density in range_density:
        graph_out_dir = Path('data_syn/density_generated_graphs')
        if not graph_out_dir.exists():
            graph_out_dir.mkdir()

        graph_out_file = Path(str(density) + '_graph.edges')
        args = ('../cmake-build-debug/bter_run',
                '-n', str(number_of_vertices),
                '-b', str(max_degree_bound),
                '-d', str(average_degree),
                '-m', str(max_cluster_coefficient),
                '-e', str(density),
                '--gpu',
                '-f', str(graph_out_dir / graph_out_file))

        print(str(density))
        popen = subprocess.Popen(args, stdout=subprocess.PIPE)
        popen.wait()
        output = popen.stdout.read()
        print(output)

    for shortest_path in range_shortest_path:
        graph_out_dir = Path('data_syn/shortest_path_generated_graphs')
        if not graph_out_dir.exists():
            graph_out_dir.mkdir()

        graph_out_file = Path(str(shortest_path) + '_graph.edges')
        args = ('../cmake-build-debug/bter_run',
                '-n', str(number_of_vertices),
                '-b', str(max_degree_bound),
                '-d', str(average_degree),
                '-m', str(max_cluster_coefficient),
                '-s', str(shortest_path),
                '--gpu',
                '-f', str(graph_out_dir / graph_out_file))

        print(str(shortest_path))
        popen = subprocess.Popen(args, stdout=subprocess.PIPE)
        popen.wait()
        output = popen.stdout.read()
        print(output)


if __name__ == '__main__':
    test_graphs_from_real_graphs()
