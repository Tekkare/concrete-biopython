from Bio.Seq import Seq as BioSeq
from Bio.Seq import MutableSeq as BioMutableSeq


_letters_to_integers_mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
_integers_to_letters_mapping = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

def _toIntegers(seq):
    string_seq = seq.__str__()

    # Verify string
    valid_chars = set(_letters_to_integers_mapping.keys())
    if not set(string_seq).issubset(valid_chars):
        raise ValueError("Invalid characters, must be A,C,G,T")
    
    # Map letters to integers
    return np.array([_letters_to_integers_mapping[char] for char in string_seq])        


def _fromIntegers(integers_array, _class):
    # Verify input array
    valid_integers = set(_integers_to_letters_mapping.keys())
    if not set(integers_array).issubset(valid_integers):
        raise ValueError("Invalid integers, must be 0,1,2,3")
    
    # Map integers to letters
    string_seq = ''.join([_integers_to_letters_mapping[integer] for integer in integers_array])

    return _class(string_seq)



class Seq(BioSeq):

    def __init__(self, data, length=None):
        super().__init__(data, length)

    def toIntegers(self):
        return _toIntegers(self)

    def fromIntegers(integers_array):
        return _fromIntegers(integers_array, Seq)


class MutableSeq(BioMutableSeq):

    def __init__(self, data, length=None):
        super().__init__(data, length)

    def toIntegers(self):
        return _toIntegers(self)

    def fromIntegers(integers_array):
        return _fromIntegers(integers_array, MutableSeq)



### tests

import unittest
import numpy as np

class MyTestCase(unittest.TestCase):

    def test_toIntegers(self):
        assert( np.all(np.equal(Seq('ACGT').toIntegers(),np.array([0,1,2,3])) ))
        assert( np.all(np.equal(Seq('TTTT').toIntegers(),np.array([3,3,3,3])) ))
        assert( np.all(np.equal(Seq('').toIntegers(),np.array([])) ))
        with self.assertRaises(ValueError):
            Seq('ACGX').toIntegers()

    def test_map_integers_to_letters(self):
        assert( Seq.fromIntegers([0,1,2,3]).__str__() == 'ACGT')
        assert( Seq.fromIntegers([3,3,3,3]).__str__() == 'TTTT')
        assert( Seq.fromIntegers([]).__str__() == '')
        with self.assertRaises(ValueError):
            Seq.fromIntegers([0,1,2,4])

if __name__ == "__main__":
    unittest.main()
