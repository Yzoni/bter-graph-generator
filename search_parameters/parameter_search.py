import scipy.optimize
import numpy as np

def evaluate_degree_distribution(alpha, beta, maxdeg, bnd, avgdeg):
    p = discrete_generalized_log_normal_probability(maxdeg, alpha, beta)
    if p[-1] > bnd:
        y1 = np.power((np.exp(1 + p[-1] - bnd)), 2) - 1
    else:
        y1 = 0

    a = np.dot(np.linspace(1, maxdeg, maxdeg), p)
    y2 = np.power(np.subtract(a, avgdeg), 2)

    y = y1 + y2

    return y


def discrete_generalized_log_normal_probability(n, alpha, beta):
    """
    Discrete generalized log-normal probability density function

    :param n:
    :param alpha:
    :param beta:
    :return:
    """
    N = np.log(np.transpose(np.linspace(1, n, n)))
    p = np.exp(-np.power(np.divide(N, alpha), beta))
    return np.divide(p, np.sum(p))


def degree_distribution_parameter_search(maxdeg, tau, avgdeg):
    handle = lambda x: evaluate_degree_distribution(x[0], x[1], maxdeg, tau, avgdeg)
    xopt = scipy.optimize.fmin(func=handle, x0=[2, 2])
    return xopt


def generate_degree_distribution(n, pdf, cutoff=0):
    """
    Create a random degree distribution from a given PDF

    :param n: Number of nodes
    :param pdf: Probability density function
    :param cutoff:
    :return:
    """
    n1 = np.sum(cutoff) # Todo
    dd1 = []


    n2 = n - n1
    tail_pdf = pdf[cutoff:] / np.sum(pdf[cutoff:])
    tail_cdf = np.cumsum(tail_pdf)
    idx2 = np.where(tail_cdf < 1)[0][-1]
    tail_cdf = np.insert(tail_cdf[:idx2], 0, 0)
    tail_cdf = np.append(tail_cdf, 1)
    coins = np.random.rand(n2, 1)

    cnts = np.histogram(coins, tail_cdf)

    idx3 = np.where(cnts[0] > 0)[0][-1]
    dd2 = cnts[0][:idx3]

    # return concatination of dd1 and dd2
    return np.append(dd1, dd2)

def optimal_xi(nd):
    max_degree = np.where(nd > 0.1)[0][-1]

    pass