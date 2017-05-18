from parameter_search import degree_distribution_parameter_search
from parameter_search import discrete_generalized_log_normal_probability
from parameter_search import generate_degree_distribution

if __name__ == '__main__':
    nnodes = 20
    maxdeg_bound = 10
    avgdeg_target = 3
    maxccd_target = 0.95
    gcc_target = 0.15
    tau = 1e-3 / nnodes

    alpha, beta = degree_distribution_parameter_search(maxdeg_bound, tau, avgdeg_target)
    pdf = discrete_generalized_log_normal_probability(nnodes, alpha, beta)
    nd = generate_degree_distribution(nnodes, pdf)
