import networkx as nx
import argparse
import glob
from pathlib import Path
import csv


def iterate_edge_lists(directory_path: Path) -> nx.Graph and Path:
    if directory_path.is_dir():
        for file_name in glob.iglob(str(directory_path) + '/**/*.edges', recursive=True):
            file_path = Path(file_name)
            with open(str(file_path), 'r') as c:
                dialect = csv.Sniffer().sniff(c.read(1024))
                csv_options = {'delimiter': dialect.delimiter, 'quotechar': dialect.quotechar}

            gx = nx.read_edgelist(str(file_path), delimiter=csv_options['delimiter'], create_using=nx.Graph(),
                                  nodetype=int)

            yield gx, file_path


def extract_largest_component(graph: nx.Graph) -> nx.Graph:
    return max(nx.connected_component_subgraphs(graph), key=len)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extraction of diversity of graph properties')
    parser.add_argument('data_input', metavar='I', type=str,
                        help='Input data directory containing edge lists')
    parser.add_argument('data_output', metavar='O', type=str,
                        help='Output directory for fixed edges lists')

    args = parser.parse_args()

    input_path = Path(args.data_input)
    out_path = Path(args.data_output)
    for edge_list, file_name in iterate_edge_lists(input_path):
        print(file_name)
        gx = extract_largest_component(edge_list)
        if not Path(args.data_output).exists():
            out_path.mkdir()
        nx.write_edgelist(gx, str(out_path / Path(str(file_name.stem) + '.' + str(file_name.suffix))),
                          delimiter=';', data=False)
