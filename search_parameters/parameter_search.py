import numpy as np
import scipy.optimize


class ParameterSearch:
    def __init__(self, number_of_nodes, max_deg_bound, avg_deg_target, max_ccd_target, gcc_target):
        self.n = number_of_nodes
        self.max_deg_bound = max_deg_bound
        self.avg_deg_target = avg_deg_target
        self.max_ccd_target = max_ccd_target
        self.gcc_target = gcc_target

        self.ccd = []
        self.nd = []

    def evaluate_degree_distribution(self, alpha, beta, bnd, avgdeg):
        p = self.discrete_generalized_log_normal_probability(alpha, beta)

        if p[-1] > bnd:
            y1 = np.power((np.exp(1 + p[-1] - bnd)), 2) - 1
        else:
            y1 = 0

        a = np.dot(np.linspace(1, self.max_deg_bound, self.max_deg_bound), p)
        y2 = np.power(np.subtract(a, avgdeg), 2)

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

    def degree_distribution_parameter_search(self, tau, avgdeg):
        handle = lambda x: self.evaluate_degree_distribution(x[0], x[1], tau, avgdeg)
        xopt = scipy.optimize.fmin(func=handle, x0=[2, 2])
        return xopt

    def generate_degree_distribution(self, pdf, cutoff=0):
        """
        Create a random degree distribution from a given PDF
    
        :param pdf: Probability density function
        :param cutoff:
        :return:
        """
        n1 = np.sum(cutoff)  # Todo
        dd1 = []

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

    def optimal_xi(self, nd):
        max_degree = np.where(nd > 0.1)[0][-1]

        pass

    def write_to_file(self, file_name):
        with open(file_name, 'w') as f:
            f.write(self.n + '\n')
            for e in self.nd:
                f.write(e + ' ')
            f.write('\n')
            for e in self.ccd:
                f.write(e + ' ')
            f.write('\n')
