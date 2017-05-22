import unittest
import numpy as np

from parameter_search import ParameterSearch


class TestParameters(unittest.TestCase):
    def test_dglnpdf(self):
        number_of_nodes = 20
        max_deg_bound = 10
        avg_deg_target = 3
        max_ccd_target = 0.95
        gcc_target = 0.15
        tau = 1e-3 / number_of_nodes
        p = ParameterSearch(number_of_nodes, max_deg_bound, avg_deg_target, max_ccd_target, gcc_target, tau)


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
                                     p.discrete_generalized_log_normal_probability(alpha, beta))


    def test_evaluate_degree_distribution(self):
        number_of_nodes = 20
        max_deg_bound = 10
        avg_deg_target = 3
        max_ccd_target = 0.95
        gcc_target = 0.15
        tau = 1e-3 / number_of_nodes

        p = ParameterSearch(number_of_nodes, max_deg_bound, avg_deg_target, max_ccd_target, gcc_target, tau)

        alpha, beta = 2, 2
        self.assertEqual(8.763245398956006, p.evaluate_degree_distribution(alpha, beta))


    def test_alpha_beta(self):
        number_of_nodes = 20
        max_deg_bound = 10
        avg_deg_target = 3
        max_ccd_target = 0.95
        gcc_target = 0.15
        tau = 1e-3 / number_of_nodes

        p = ParameterSearch(number_of_nodes, max_deg_bound, avg_deg_target, max_ccd_target, gcc_target, tau)

        alpha, beta = p.degree_distribution_parameter_search()
        self.assertAlmostEqual(1.717337, np.round(alpha, 6))
        self.assertAlmostEqual(7.945021, np.round(beta, 6))

    def test_find_xi(self):
        pass


if __name__ == '__main__':
    unittest.main()
