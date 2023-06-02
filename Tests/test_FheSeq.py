import sys, os
sys.path.append(os.getcwd())

from BioConcrete.FheSeq import *
from BioConcrete.utils import *

import unittest
import numpy as np
from Bio.Seq import Seq, MutableSeq
from concrete import fhe

class BioConcreteCircuit:
   
    def __init__(self, seq_length):
        self.seq_length=seq_length
        self.inputset=[
            (np.random.randint(0, 4, size=(seq_length,)),
            np.random.randint(0, 4, size=(seq_length,)))
            for _ in range(100)
        ]
        self.circuit = None

    def set(self, circuitFunction, verbose=False):
        compiler = fhe.Compiler(lambda data1,data2: circuitFunction(data1, data2), {"data1": "encrypted", "data2": "encrypted"})
        self.circuit = compiler.compile(
            inputset=self.inputset,
            configuration=fhe.Configuration(
                enable_unsafe_features=True,
                use_insecure_key_cache=True,
                insecure_key_cache_location=".keys",
                #dataflow_parallelize=True,
            ),
            verbose=verbose,
        )
    
    def run(self, seq1, seq2, simulate=False, integers_output=False):
        if not self.circuit: raise Error('circuit was not set')
        assert len(seq1) == self.seq_length, f"Sequence 1 length is not correct, should be {self.seq_length} characters"
        assert len(seq2) == self.seq_length, f"Sequence 2 length is not correct, should be {self.seq_length} characters"

        # convert letters to integers
        integers1 = seqToIntegers(seq1)
        integers2 = seqToIntegers(seq2)
        output_seq = self.circuit.simulate(integers1, integers2) if simulate else self.circuit.encrypt_run_decrypt(integers1, integers2)

        if not integers_output:
            # convert back integers to letters    
            return seqFromIntegers(output_seq)
        else:
            return output_seq


class TestFheSeq(unittest.TestCase):

    def test_operands(self):
        circuit = BioConcreteCircuit(4)
        seq1 = Seq('ACGT')    
        seq2 = Seq('CGTA')        
        seq3 = Seq('ACGA')        

        # == operands
        circuit.set(lambda x,y: FheSeq(x)==FheSeq(y) )
        assert( np.all(circuit.run(seq1, seq1, True, True)) )        
        assert( not np.all(circuit.run(seq1, seq3, True, True)) )    

        # < operands
        circuit.set(lambda x,y: FheSeq(x)<FheSeq(y) )
        assert( np.all(circuit.run(seq1, seq2, True, True) == np.array([True, True, True, False])) )  

        # len operand
        circuit.set( lambda x,y: fhe.ones(1) if len(FheSeq(x))==4 else fhe.zeros(1) )
        assert( np.all(circuit.run(seq1, seq2, True, True) == np.array([1])) )

        # single getitem
        def getitem(x,y):
            out=fhe.zeros(4)
            eseq = FheSeq(x)
            for i in range(4):
                out[i]=eseq[i]
            return out
        circuit.set(getitem)
        assert( circuit.run(seq1, seq1, True, False) == seq1 )

        # multiple getitem
        def getitems(x,y):
            eseq = FheSeq(x)
            return eseq[0:2].toArray()
        circuit.set(getitems)
        assert( circuit.run(seq1, seq1, True, False) == seq1[0:2] )

        # # add sequences (concat them)
        # circuit.set( lambda x,y: FheSeq(x)+FheSeq(y) )
        # assert( circuit.run(seq1, seq2, True, False) == (seq1+seq2) )

        # # add array to the left (concat them)
        # circuit.set( lambda x,y: fhe.zeros(2)+FheSeq(x) )
        # assert( circuit.run(seq1, seq2, True, False) == ('AA'+seq1) )        

    def test_complement(self):
        circuit = BioConcreteCircuit(4)
        seq1 = Seq('ACGT')        
        circuit.set(lambda x,y:FheSeq(x).complement().toArray())
        assert( circuit.run(seq1, seq1, True) == seq1.complement())

    def test_reverse_complement(self):
        circuit = BioConcreteCircuit(4)
        seq1 = Seq('ACGT')        
        circuit.set(lambda x,y:FheSeq(x).reverse_complement().toArray())
        assert( circuit.run(seq1, seq1, True) == seq1.reverse_complement())



if __name__ == "__main__":
    #unittest.main()

    suite = unittest.TestLoader().loadTestsFromName('test_FheSeq.TestFheSeq.test_operands')
    unittest.TextTestRunner(verbosity=2).run(suite)    
