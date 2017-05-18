import unittest
from parameter_search import degree_distribution_parameter_search
from parameter_search import evaluate_degree_distribution
from parameter_search import discrete_generalized_log_normal_probability
from parameter_search import generate_degree_distribution
import numpy as np


class TestParameters(unittest.TestCase):
    def test_dglnpdf(self):
        maxdeg_bound = 10
        alpha = 2
        beta = 2
        np.testing.assert_array_almost_equal(np.transpose([0.181540623715284,
                                                           0.160993855120665,
                                                           0.134255165662288,
                                                           0.112283445406358,
                                                           0.095002846307965,
                                                           0.081359609485210,
                                                           0.070445395213683,
                                                           0.061587576466916,
                                                           0.054300340083636,
                                                           0.048231142537995]),
                                             discrete_generalized_log_normal_probability(maxdeg_bound, alpha, beta))

    def test_dglnobjfunc(self):
        self.assertEqual(8.763245398956006, evaluate_degree_distribution(2, 2, 10, 5e-5, 3))

    def test_alpha_beta(self):
        nnodes = 20
        maxdeg_bound = 10
        avgdeg_target = 3
        tau = 1e-3 / nnodes
        alpha, beta = degree_distribution_parameter_search(maxdeg_bound, tau, avgdeg_target)
        self.assertAlmostEqual(1.717337, np.round(alpha, 6))
        self.assertAlmostEqual(7.945021, np.round(beta, 6))

    def test_generate_degree_distribution(self):
        nnodes = 20
        maxdeg_bound = 10
        alpha = 2
        beta = 2
        pdf = discrete_generalized_log_normal_probability(maxdeg_bound, alpha, beta)
        print(generate_degree_distribution(nnodes, pdf))

    def test_find_xi(self):
        pass


if __name__ == '__main__':
    unittest.main()
