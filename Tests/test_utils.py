import sys, os
sys.path.append(os.getcwd())

from BioConcrete.utils import *

import unittest
import numpy as np
from Bio.Seq import Seq, MutableSeq

class TestUtils(unittest.TestCase):

    def test_seqToIntegers(self):
        assert( np.all(np.equal(seqToIntegers(Seq('ACGT')),np.array([0,1,2,3])) ))
        assert( np.all(np.equal(seqToIntegers(Seq('TTTT')),np.array([3,3,3,3])) ))
        assert( np.all(np.equal(seqToIntegers(Seq('')),np.array([])) ))
        with self.assertRaises(ValueError):
            seqToIntegers(Seq('ACGX'))

    def test_seqFromIntegers(self):
        assert( seqFromIntegers([0,1,2,3]).__str__() == 'ACGT')
        assert( seqFromIntegers([3,3,3,3]).__str__() == 'TTTT')
        assert( seqFromIntegers([]).__str__() == '')
        with self.assertRaises(ValueError):
            seqFromIntegers([0,1,2,4])

if __name__ == "__main__":
    unittest.main()
