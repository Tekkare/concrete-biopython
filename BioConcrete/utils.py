import numpy as np
from Bio.Seq import Seq, MutableSeq

_letters_to_integers_mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
_integers_to_letters_mapping = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}


def seqToIntegers(seq):
    string_seq = seq.__str__()

    # Verify string
    valid_chars = set(_letters_to_integers_mapping.keys())
    if not set(string_seq).issubset(valid_chars):
        raise ValueError("Invalid characters, must be A,C,G,T")
    
    # Map letters to integers
    return np.array([_letters_to_integers_mapping[char] for char in string_seq])        


def _seqFromIntegers(integers_array, _class):
    # Verify input array
    valid_integers = set(_integers_to_letters_mapping.keys())
    if not set(integers_array).issubset(valid_integers):
        raise ValueError("Invalid integers, must be 0,1,2,3")
    
    # Map integers to letters
    string_seq = ''.join([_integers_to_letters_mapping[integer] for integer in integers_array])

    return _class(string_seq)


def seqFromIntegers(integers_array):
    return _seqFromIntegers(integers_array, Seq)

def mutableSeqFromIntegers(integers_array):
    return _seqFromIntegers(integers_array, MutableSeq)