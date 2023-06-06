import sys, os
sys.path.append(os.getcwd())

from concreteBiopython.FheSeq import FheSeq, FheMutableSeq
from concreteBiopython.SeqWrapper import SeqWrapper

import unittest
import numpy as np
from Bio.Seq import Seq, MutableSeq
from concrete import fhe


class BioConcreteCircuit:
   
    def __init__(self, seq_length, simulate=False):
        self.seq_length=seq_length
        self.inputset=[
            (np.random.randint(0, len(SeqWrapper.LETTERS), size=(seq_length,)),
            np.random.randint(0, len(SeqWrapper.LETTERS), size=(seq_length,)))
            for _ in range(100)
        ]
        self.circuit = None
        self.simulate = simulate

    def set(self, circuitFunction, integers_output=False, verbose=False):
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
        self._integers_output = integers_output
    
    def run(self, seq1, seq2):
        if not self.circuit: raise Error('circuit was not set')
        assert len(seq1) == self.seq_length, f"Sequence 1 length is not correct, should be {self.seq_length} characters"
        assert len(seq2) == self.seq_length, f"Sequence 2 length is not correct, should be {self.seq_length} characters"

        # convert letters to integers
        integers1 = SeqWrapper(seq1).toIntegers()
        integers2 = SeqWrapper(seq2).toIntegers()
        int_output = self.circuit.simulate(integers1, integers2) if self.simulate else self.circuit.encrypt_run_decrypt(integers1, integers2)

        if not self._integers_output:
            # convert back integers to letters    
            return SeqWrapper(int_output).toSeq()
        else:
            return int_output


SIMULATE=False

class TestFheSeq(unittest.TestCase):

    def test_operands(self):
        circuit = BioConcreteCircuit(5, SIMULATE)
        seq1 = Seq('ACGTU')    
        seq2 = Seq('CGTUA')        
        seq3 = Seq('ACGTA')        

        # == operands
        circuit.set(lambda x,y: FheSeq(x)==FheSeq(y) , True)
        assert( circuit.run(seq1, seq1) )
        assert( not circuit.run(seq1, seq3))

        # # < operands
        # circuit.set(lambda x,y: FheSeq(x)<FheSeq(y) , True)
        # assert( np.all(circuit.run(seq1, seq2) == np.array([True, True, True, True, False])) )  

        # len operand
        circuit.set( lambda x,y: fhe.ones(1) if len(FheSeq(x))==5 else fhe.zeros(1), True )
        assert( np.all(circuit.run(seq1, seq2) == np.array([1])) )

        # single getitem
        def getitem(x,y):
            out=fhe.zeros(5)
            eseq = FheSeq(x)
            for i in range(5):
                out[i]=eseq[i]
            return out
        circuit.set(getitem)
        assert( circuit.run(seq1, seq1) == seq1 )

        # multiple getitem
        def getitems(x,y):
            eseq = FheSeq(x)
            return eseq[0:2].toArray()
        circuit.set(getitems)
        assert( circuit.run(seq1, seq1) == seq1[0:2] )

        # add sequences (concat them)
        circuit.set( lambda x,y: (FheSeq(x)+FheSeq(y)).toArray() )
        assert( circuit.run(seq1, seq2) == (seq1+seq2) )

    def test_startswith(self):
        circuit = BioConcreteCircuit(4, SIMULATE)
        seq1 = Seq('ACGT')        
        seq2 = Seq('CGTA')
        circuit.set(lambda x,y:FheSeq(x).startswith(FheSeq(y)[0:2]), True)
        assert( circuit.run(seq1, seq1) )
        assert( not circuit.run(seq1, seq2) )

        # also test with naked array:
        circuit.set(lambda x,y:FheSeq(x).startswith(FheSeq(y)[0:2].toArray()), True)
        assert( circuit.run(seq1, seq1) )

        # also test with start or end
        circuit.set(lambda x,y:FheSeq(x).startswith(FheSeq(y)[0:2],1), True)
        assert( circuit.run(seq1, seq2) )
        circuit.set(lambda x,y:FheSeq(x).startswith(FheSeq(y)[0:2],None, 3), True)
        assert( circuit.run(seq1, seq2) )   

    def test_endswith(self):
        circuit = BioConcreteCircuit(4, SIMULATE)
        seq1 = Seq('ACGT')        
        seq2 = Seq('CGTA')
        circuit.set(lambda x,y:FheSeq(x).endswith(FheSeq(y)[1:3]), True)
        assert( not circuit.run(seq1, seq1) )
        assert( circuit.run(seq1, seq2) )

    def test_complement(self):
        circuit = BioConcreteCircuit(4, SIMULATE)
        seq1 = Seq('ACGT')        
        circuit.set(lambda x,y:FheSeq(x).complement().toArray())
        assert( circuit.run(seq1, seq1) == seq1.complement())

        with self.assertRaises(TypeError): # verify we cannot change a immutable seq in place
            circuit.set(lambda x,y:FheSeq(x).complement(True).toArray())
            circuit.run(seq1, seq1)

        # test changing mutable seq in place
        seq2 = MutableSeq('ACGT')
        circuit.set(lambda x,y:FheMutableSeq(x).complement(True).toArray())
        assert( circuit.run(seq2, seq2) == seq2.complement(True))

    def test_complement_rna(self):
        circuit = BioConcreteCircuit(4, SIMULATE)
        seq1 = Seq('ACGU')        
        circuit.set(lambda x,y:FheSeq(x).complement_rna().toArray())
        assert( circuit.run(seq1, seq1) == seq1.complement_rna())

        with self.assertRaises(TypeError): # verify we cannot change a immutable seq in place
            circuit.set(lambda x,y:FheSeq(x).complement_rna(True).toArray())
            circuit.run(seq1, seq1)

        # test changing mutable seq in place
        seq2 = MutableSeq('ACGU')
        circuit.set(lambda x,y:FheMutableSeq(x).complement(True).toArray())
        assert( circuit.run(seq2, seq2) == seq2.complement(True))            

    def test_reverse_complement(self):
        circuit = BioConcreteCircuit(4, SIMULATE)
        seq1 = Seq('ACGT')        
        circuit.set(lambda x,y:FheSeq(x).reverse_complement().toArray())
        assert( circuit.run(seq1, seq1) == seq1.reverse_complement())

        with self.assertRaises(TypeError): # verify we cannot change a immutable seq in place
            circuit.set(lambda x,y:FheSeq(x).reverse_complement(True).toArray())
            circuit.run(seq1, seq1)

    def test_reverse_complement_rna(self):
        circuit = BioConcreteCircuit(4, SIMULATE)
        seq1 = Seq('ACGU')        
        circuit.set(lambda x,y:FheSeq(x).reverse_complement_rna().toArray())
        assert( circuit.run(seq1, seq1) == seq1.reverse_complement_rna())

        with self.assertRaises(TypeError): # verify we cannot change a immutable seq in place
            circuit.set(lambda x,y:FheSeq(x).reverse_complement_rna(True).toArray())
            circuit.run(seq1, seq1)

    def test_transcribe(self):
        circuit = BioConcreteCircuit(4, SIMULATE)
        seq1 = Seq('ACGT')        
        circuit.set(lambda x,y:FheSeq(x).transcribe().toArray())
        assert( circuit.run(seq1, seq1) == seq1.transcribe())

        with self.assertRaises(TypeError): # verify we cannot change a immutable seq in place
            circuit.set(lambda x,y:FheSeq(x).transcribe(True).toArray())
            circuit.run(seq1, seq1)

    def test_back_transcribe(self):
        circuit = BioConcreteCircuit(4, SIMULATE)
        seq1 = Seq('ACGU')        
        circuit.set(lambda x,y:FheSeq(x).back_transcribe().toArray())
        assert( circuit.run(seq1, seq1) == seq1.back_transcribe())

        with self.assertRaises(TypeError): # verify we cannot change a immutable seq in place
            circuit.set(lambda x,y:FheSeq(x).back_transcribe(True).toArray())
            circuit.run(seq1, seq1)   

    def test_translate(self):
        circuit = BioConcreteCircuit(39, SIMULATE)
        seq1 = Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')        
        circuit.set(lambda x,y:FheSeq(x).translate().toArray())
        assert( circuit.run(seq1, seq1) == seq1.translate())

    # def test_join(self):
    #     seq1 = Seq('ACG')
    #     seq2 = Seq('TTT')
    #     circuit = BioConcreteCircuit(3, SIMULATE)
    #     circuit.set(lambda x,y:FheSeq(x).join(FheSeq(y)).toArray())
    #     assert( circuit.run(seq2, seq1) == seq2.join(seq1))

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--sim', action='store_true', help='simulate encryption')
args = parser.parse_args()


#SIMULATE = args.sim
SIMULATE = True
#unittest.main()

suite = unittest.TestLoader().loadTestsFromName('test_FheSeq.TestFheSeq.test_join')
unittest.TextTestRunner(verbosity=2).run(suite)    
