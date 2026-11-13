import unittest
from export_t_sequences import scaled_ranks
class RankTests(unittest.TestCase):
    def test_ties(self):
        self.assertEqual(scaled_ranks([1, 1, 2, 3]), {1:0, 2:.6, 3:1})
    def test_order(self):
        self.assertEqual(scaled_ranks([3, 1, 2]), {1:0, 2:.5, 3:1})
    def test_constant_rejected(self):
        with self.assertRaises(ValueError): scaled_ranks([1, 1])
if __name__ == '__main__': unittest.main()
