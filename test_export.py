import unittest

from export_t_sequences import scaled_ranks


class ScaledRanksTests(unittest.TestCase):
    def test_ties_share_average_rank(self):
        self.assertEqual(scaled_ranks([1, 1, 2, 3]), {1: 0, 2: 0.6, 3: 1})

    def test_ranks_follow_value_order(self):
        self.assertEqual(scaled_ranks([3, 1, 2]), {1: 0, 2: 0.5, 3: 1})

    def test_constant_target_rejected(self):
        with self.assertRaises(ValueError):
            scaled_ranks([1, 1])


if __name__ == '__main__':
    unittest.main()
