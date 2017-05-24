from parameters.search import ParameterSearch

if __name__ == '__main__':
    number_of_nodes = 20
    max_deg_bound = 10
    avg_deg_target = 3
    max_ccd_target = 0.95
    gcc_target = 0.15
    tau = 1e-3 / number_of_nodes

    print(number_of_nodes)

    p = ParameterSearch(number_of_nodes, max_deg_bound, avg_deg_target, max_ccd_target, gcc_target, tau)
    nd = p.run_nd(number_of_nodes, max_deg_bound, avg_deg_target, max_ccd_target, gcc_target, tau)
    print('nd: ' + str(nd))
    # print('ccd: ' + str(ccd))

    # ParameterSearch.write_to_file('parameters.txt', nd, ccd)
