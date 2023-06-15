import sys, os

import unittest
import numpy as np
from concrete import fhe
from Bio.Seq import Seq, MutableSeq

sys.path.append(os.getcwd())

from concrete_biopython.FheSeq import FheSeq, FheMutableSeq
from concrete_biopython.SeqWrapper import SeqWrapper


class BioConcreteCircuit:
   
    """
    Circuit factory class for testing FheSeq on 2 sequences input
    """
    def __init__(self, seq_length, simulate=False):
        self.seq_length=seq_length
        self.inputset=[
            (np.random.randint(0, SeqWrapper.maxInteger()+1, size=(seq_length,)),
            np.random.randint(0, SeqWrapper.maxInteger()+1, size=(seq_length,)))
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

        ## == operand
        circuit.set(lambda x,y: FheSeq(x)==FheSeq(y) , True)
        assert( circuit.run(seq1, seq1) )
        assert( not circuit.run(seq1, seq3))

        ## >= operand
        circuit.set(lambda x,y: FheSeq(x)>=FheSeq(y) , True)

        # close sequences:
        assert( circuit.run(seq1, seq3) )        
        assert( not circuit.run(seq3, seq1) )

        assert( circuit.run(seq2, seq3) )        
        assert( not circuit.run(seq3, seq2) )

        # identical sequences
        assert( circuit.run(seq1, seq1) )

        # sequences with different sizes:
        seq4 = Seq('AAAAA')
        circuit.set(lambda x,y: FheSeq(x)>=FheSeq(y)[0:4] , True)        
        assert( circuit.run(seq4, seq4) )
        circuit.set(lambda x,y: FheSeq(x)[0:4]>=FheSeq(y) , True)        
        assert( not circuit.run(seq4, seq4) )


        ## <= operand
        circuit.set(lambda x,y: FheSeq(x)<=FheSeq(y) , True)

        # close sequences:
        assert( circuit.run(seq3, seq1) )        
        assert( not circuit.run(seq1, seq3) )

        assert( circuit.run(seq3, seq2) )        
        assert( not circuit.run(seq2, seq3) )

        # identical sequences
        assert( circuit.run(seq1, seq1) )

        # sequences with different sizes:
        seq4 = Seq('AAAAA')
        circuit.set(lambda x,y: FheSeq(x)<=FheSeq(y)[0:4] , True)        
        assert( not circuit.run(seq4, seq4) )
        circuit.set(lambda x,y: FheSeq(x)[0:4]<=FheSeq(y) , True)        
        assert( circuit.run(seq4, seq4) )


        ## > operand
        circuit.set(lambda x,y: FheSeq(x)>FheSeq(y) , True)

        # close sequences:
        assert( circuit.run(seq1, seq3) )        
        # identical sequences
        assert( not circuit.run(seq1, seq1) )
        # sequences with different sizes:
        circuit.set(lambda x,y: FheSeq(x)>FheSeq(y)[0:4] , True)        
        assert( circuit.run(seq4, seq4) )

        ## < operand
        circuit.set(lambda x,y: FheSeq(x)<FheSeq(y) , True)

        # close sequences:
        assert( not circuit.run(seq1, seq3) )        
        # identical sequences
        assert( not circuit.run(seq1, seq1) )
        # sequences with different sizes:
        circuit.set(lambda x,y: FheSeq(x)[0:4]<FheSeq(y) , True)        
        assert( circuit.run(seq4, seq4) )        

        ## >= operand
        circuit.set(lambda x,y: FheSeq(x)>=FheSeq(y) , True)

        # close sequences:
        assert( circuit.run(seq1, seq3) )        

        ## len operand
        circuit.set( lambda x,y: fhe.ones(1) if len(FheSeq(x))==5 else fhe.zeros(1), True )
        assert( np.all(circuit.run(seq1, seq2) == np.array([1])) )

        ## single getitem
        def getitem(x,y):
            out=fhe.zeros(5)
            eseq = FheSeq(x)
            for i in range(5):
                out[i]=eseq[i]
            return out
        circuit.set(getitem)
        assert( circuit.run(seq1, seq1) == seq1 )

        ## multiple getitem
        def getitems(x,y):
            eseq = FheSeq(x)
            return eseq[0:2].toArray()
        circuit.set(getitems)
        assert( circuit.run(seq1, seq1) == seq1[0:2] )

        ## multiple getitem with copy
        def getitemscp(x,y):
            eseq = FheMutableSeq(x)
            eseq2= eseq[0:2]
            eseq2[0] = eseq2[1]
            return eseq[0:2].toArray()
        circuit.set(getitemscp)
        assert( circuit.run(seq1, seq1) == seq1[0:2] )        

        # add sequences (concat them)
        circuit.set( lambda x,y: (FheSeq(x)+FheSeq(y)).toArray() )
        assert( circuit.run(seq1, seq2) == (seq1+seq2) )

        ## >= operand for 1 letter
        circuit2 = BioConcreteCircuit(1, SIMULATE)
        seq1_2 = Seq('C')
        seq2_2 = Seq('A')

        circuit2.set(lambda x,y: FheSeq(x)>=FheSeq(y) , True)
        print( 'RES :',circuit2.run(seq1_2, seq2_2))
        assert( circuit2.run(seq1_2, seq2_2) )
        assert( circuit2.run(seq1_2, seq1_2) )
        assert( not circuit2.run(seq2_2, seq1_2) )

    def test_iter(self):
        seq1 = Seq('ACGT')
        seq2 = Seq('AACGT')
        circuit = BioConcreteCircuit(4, SIMULATE)
       
        def iter(x,y):
            seq=FheSeq(x)
            seq2=FheMutableSeq(seq[0])
            for c in seq:
                seq2.append(c)
            return seq2.toArray()
        circuit.set( iter ) 
        assert( circuit.run(seq1, seq1) == seq2) 

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
        circuit.set(lambda x,y:FheSeq(x).startswith(FheSeq(y)[3],None, 3), True)
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

    def test_join(self):
        seq1 = Seq('ACG')
        seq2 = Seq('TTT')
        circuit = BioConcreteCircuit(3, SIMULATE)
       
        # test joining a FheSeq
        circuit.set(lambda x,y:FheSeq(x).join(FheSeq(y)).toArray())
        assert( circuit.run(seq2, seq1) == seq2.join(seq1))
        
        # test joining an array
        circuit.set(lambda x,y:FheSeq(x).join(y).toArray())
        assert( circuit.run(seq2, seq1) == seq2.join(seq1))
        
        # test joining a sequence of FheSeq and arrays
        circuit.set(lambda x,y:FheSeq(x).join([FheSeq(x),x,y]).toArray())
        assert( circuit.run(seq2, seq1) == seq2.join([seq2,seq2,seq1]))   

    def test_setitem(self):
        seq1 = Seq('XXXXXXXX')
        seq2 = Seq('ABCDEFGH')
        circuit = BioConcreteCircuit(8, SIMULATE)
       
        def setitem(x,y):
            seq=FheMutableSeq(x)
            seq[0]=y[0]
            seq[2:4]=y[2:4]
            seq[5:]=FheSeq(y[5:])
            return seq.toArray()
        circuit.set( setitem )
        assert( circuit.run(seq1, seq2) == Seq('AXCDXFGH'))

    def test_delitem(self):
        seq1 = Seq('ABCDEFGH')
        seq2 = MutableSeq('ABCDEFGH')
        circuit = BioConcreteCircuit(8, SIMULATE)
       
        def delitem(x,y):
            seq=FheMutableSeq(x)
            del seq[0]
            del seq[2]
            del seq[2:4]
            return seq.toArray()
        circuit.set( delitem )
        del seq2[0]
        del seq2[2]
        del seq2[2:4]
        assert( circuit.run(seq1, seq1) == seq2)

    def test_append(self):
        seq1 = Seq('AC')
        seq1_bis = MutableSeq('AC')
        seq2 = MutableSeq('GT')
        circuit = BioConcreteCircuit(2, SIMULATE)
       
        def append(x,y):
            seq=FheMutableSeq(x)
            seq2=FheMutableSeq(y)
            seq.append(seq2[0])
            return seq.toArray()
        circuit.set( append ) 
        seq1_bis.append(seq2[0])
        assert( circuit.run(seq1, seq2) == seq1_bis)  

    def test_insert(self):
        seq1 = Seq('AC')
        seq1_bis = MutableSeq('AC')
        seq2 = MutableSeq('GT')
        circuit = BioConcreteCircuit(2, SIMULATE)
       
        def insert(x,y):
            seq=FheMutableSeq(x)
            seq2=FheMutableSeq(y)
            seq.insert(1,seq2[0])
            return seq.toArray()
        circuit.set( insert ) 
        seq1_bis.insert(1,seq2[0])
        assert( circuit.run(seq1, seq2) == seq1_bis)  

    def test_pop(self):
        seq1 = Seq('ACGT')
        seq1_bis = MutableSeq('ACGT')
        circuit = BioConcreteCircuit(4, SIMULATE)
       
        def pop(x,y):
            seq=FheMutableSeq(x)
            seq.pop()
            seq.pop(1)
            return seq.toArray()
        circuit.set( pop ) 
        seq1_bis.pop()
        seq1_bis.pop(1)
        assert( circuit.run(seq1, seq1) == seq1_bis)  

    def test_extend(self):
        seq1 = Seq('AC')
        seq1_bis = MutableSeq('AC')
        seq2 = MutableSeq('GT')
        circuit = BioConcreteCircuit(2, SIMULATE)
       
        def extend(x,y):
            seq=FheMutableSeq(x)
            seq.extend(y)
            seq.extend(FheSeq(y))
            return seq.toArray()
        circuit.set( extend ) 
        seq1_bis.extend(seq2)
        seq1_bis.extend(seq2)
        assert( circuit.run(seq1, seq2) == seq1_bis)   


SIMULATE = True

unittest.main()

# suite = unittest.TestLoader().loadTestsFromName('test_FheSeq.TestFheSeq.test_operands')
# unittest.TextTestRunner(verbosity=1).run(suite)
