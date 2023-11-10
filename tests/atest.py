import numpy as np
from concrete import fhe
from Bio.Seq import Seq, MutableSeq

import sys, os, time
sys.path.append(os.getcwd())
from concrete_biopython.FheSeq import FheSeq, FheMutableSeq
from concrete_biopython.SeqWrapper import SeqWrapper
from concrete_biopython.BioCircuit import BioCircuit



def process_sequence(seq1, seq2):
    """
    Check wether two elements are identical
    This is a dummy function for testing, we will redefine it for each application later
    """
    return seq1[0]==seq2[1]


def compute_fhe_output(seq_list, process_seq, description, res=False):
    """
    Pack up the creation and running of a circuit,
    and display the result along with the computation times
    """
    # First print the description and a waiting message
    print(description)
    print("Computing ...", end="", flush=True)
    print("\r", end="")

    # prepare configuration
    configuration=fhe.Configuration(
        enable_unsafe_features=True,
        use_insecure_key_cache=True,
        insecure_key_cache_location=".keys",
        dataflow_parallelize=False,
    )  

    
    # Create a BioCircuit wrapped circuit
    circuit = BioCircuit(
        function=process_seq,
        len_seqs=[len(seq) for seq in seq_list],
        configuration=configuration,
    )

    output = circuit.encrypt_run_decrypt_as_integers(*seq_list)
    
    print('==> Result :', output, '\n')

    # In case we need the result
    if res:
        return output

# Testing the functions
seq1 = Seq('ACCAGGTAC')
seq2 = MutableSeq('CGTTAGC')

# Test process_sequence
OUTPUT_PYTHON = process_sequence(seq1, seq2)
print("Python: Are the two sequences equal ?: ", OUTPUT_PYTHON)    

# Compute our simple process_sequence function in FHE on our test sequence:
compute_fhe_output([seq1, seq2], process_sequence, description='FHE: Are the first two sequences equal ?:')
