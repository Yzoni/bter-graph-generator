import numpy as np
import scipy.optimize


class ParameterSearch:
    def __init__(self, number_of_nodes, max_deg_bound, avg_deg_target, max_ccd_target, gcc_target, tau):
        self.n = number_of_nodes
        self.max_deg_bound = max_deg_bound
        self.avg_deg_target = avg_deg_target
        self.max_ccd_target = max_ccd_target
        self.gcc_target = gcc_target
        self.tau = tau

    def evaluate_degree_distribution(self, alpha, beta):
        p = self.discrete_generalized_log_normal_probability(alpha, beta)

        if p[-1] > self.tau:
            y1 = np.power((np.exp(1 + p[-1] - self.tau)), 2) - 1
        else:
            y1 = 0

        a = np.dot(np.arange(1, self.max_deg_bound + 1, 1), p)
        y2 = np.power(np.subtract(a, self.avg_deg_target), 2)

        y = y1 + y2

        return y

    def discrete_generalized_log_normal_probability(self, alpha, beta):
        """
        Discrete generalized log-normal probability density function
    
        :param alpha:
        :param beta:
        :return:
        """
        N = np.log(np.transpose(np.linspace(1, self.max_deg_bound, self.max_deg_bound)))
        p = np.exp(-np.power(np.divide(N, alpha), beta))
        return np.divide(p, np.sum(p))

    def degree_distribution_parameter_search(self):
        handle = lambda x: self.evaluate_degree_distribution(x[0], x[1])
        xopt = scipy.optimize.fmin(func=handle, x0=[2, 2])
        return xopt

    def generate_degree_distribution(self, pdf, cutoff=0):
        """
        Create a random degree distribution from a given PDF
    
        :param pdf: Probability density function
        :param cutoff:
        :return:
        """

        dd1 = np.round(self.n * pdf[:cutoff])
        n1 = np.sum(cutoff)  # Todo

        n2 = self.n - n1
        tail_pdf = pdf[cutoff:] / np.sum(pdf[cutoff:])
        tail_cdf = np.cumsum(tail_pdf)
        idx2 = np.where(tail_cdf < 1)[0][-1]
        tail_cdf = np.insert(tail_cdf[:idx2], 0, 0)
        tail_cdf = np.append(tail_cdf, 1)
        coins = np.random.rand(n2, 1)

        cnts = np.histogram(coins, tail_cdf)

        idx3 = np.where(cnts[0] > 0)[0][-1]
        dd2 = cnts[0][:idx3]

        return np.append(dd1, dd2)

    def optimal_xi(self, nd, pdf):
        self.evaluate_degree_distribution(self.n, pdf)

        handle = lambda x: self._compute_objective(nd, x)
        xopt = scipy.optimize.fmin(func=handle, x0=[0.5])

        return xopt[0], np.where(nd > 0)[0][-1]

    def _compute_objective(self, nd, xi):
        maxd = nd.size
        ccd_mean = self.max_ccd_target * np.exp(- np.transpose(np.arange(0, maxd, 1) * xi))

        n_wedges = np.transpose(nd) * np.arange(1, maxd + 1, 1) * ((np.arange(1, maxd + 1, 1) - 1) / 2)
        gcc_xi = np.dot(n_wedges, ccd_mean) / np.sum(n_wedges)
        y = np.abs(self.gcc_target - gcc_xi)
        return y

    def target_clustering_coefficient_per_degree(self, xi, max_degree):
        return self.max_ccd_target * np.exp(-np.arange(0, max_degree + 1, 1) * xi)

    def write_to_file(self, file_name, nd, ccd):
        with open(file_name, 'w') as f:
            f.write(str(self.n) + '\n')
            for e in nd:
                f.write(str(e) + ' ')
            f.write('\n')
            for e in ccd:
                f.write(str(e) + ' ')
            f.write('\n')
