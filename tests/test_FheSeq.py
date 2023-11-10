import sys, os

import unittest
import numpy as np
from concrete import fhe
from Bio.Seq import Seq, MutableSeq

sys.path.append(os.getcwd())

from concrete_biopython.FheSeq import FheSeqMaker, Alphabets, FheSeq, FheMutableSeq
from concrete_biopython.BioCircuit import BioCircuit

FHE_SEQ_MAKER = FheSeqMaker(Alphabets.PROTEINS)

def make_circuit(function, seq_length, seq_output=False):

    # prepare configuration
    configuration=fhe.Configuration(
        enable_unsafe_features=True,
        use_insecure_key_cache=True,
        insecure_key_cache_location=".keys",
        dataflow_parallelize=False,
    )  

    # Create a BioCircuit wrapped circuit
    circuit = BioCircuit(
        function=function,
        len_seqs=[seq_length, seq_length],
        fhe_seq_maker = FHE_SEQ_MAKER,
        configuration=configuration,
        seq_output=seq_output,
        show_timing=False,
    )

    return circuit



class TestFheSeq(unittest.TestCase):

    def test_operands(self):        
        seq1 = Seq('ACGTU')    
        seq2 = Seq('CGTUA')        
        seq3 = Seq('ACGTA')        

        ## == operand
        circuit = make_circuit(lambda x,y:  x==y , 5)
        assert( circuit.encrypt_run_decrypt(seq1, seq1) )
        assert( not circuit.encrypt_run_decrypt(seq1, seq3))

        ## >= operand
        circuit = make_circuit(lambda x,y:  x>=y , 5)

        # close sequences:
        assert( circuit.encrypt_run_decrypt(seq1, seq3) )        
        assert( not circuit.encrypt_run_decrypt(seq3, seq1) )

        assert( circuit.encrypt_run_decrypt(seq2, seq3) )        
        assert( not circuit.encrypt_run_decrypt(seq3, seq2) )

        # identical sequences
        assert( circuit.encrypt_run_decrypt(seq1, seq1) )

        # sequences with different sizes:
        seq4 = Seq('AAAAA')
        circuit = make_circuit(lambda x,y:  x>=y[0:4] , 5)        
        assert( circuit.encrypt_run_decrypt(seq4, seq4) )
        circuit = make_circuit(lambda x,y:  x[0:4]>=y , 5)        
        assert( not circuit.encrypt_run_decrypt(seq4, seq4) )


        ## <= operand
        circuit = make_circuit(lambda x,y:  x<=y , 5)

        # close sequences:
        assert( circuit.encrypt_run_decrypt(seq3, seq1) )        
        assert( not circuit.encrypt_run_decrypt(seq1, seq3) )

        assert( circuit.encrypt_run_decrypt(seq3, seq2) )        
        assert( not circuit.encrypt_run_decrypt(seq2, seq3) )

        # identical sequences
        assert( circuit.encrypt_run_decrypt(seq1, seq1) )

        # sequences with different sizes:
        seq4 = Seq('AAAAA')
        circuit = make_circuit(lambda x,y:  x<=y[0:4] , 5)        
        assert( not circuit.encrypt_run_decrypt(seq4, seq4) )
        circuit = make_circuit(lambda x,y:  x[0:4]<=y , 5)        
        assert( circuit.encrypt_run_decrypt(seq4, seq4) )


        ## > operand
        circuit = make_circuit(lambda x,y:  x>y , 5)

        # close sequences:
        assert( circuit.encrypt_run_decrypt(seq1, seq3) )        
        # identical sequences
        assert( not circuit.encrypt_run_decrypt(seq1, seq1) )
        # sequences with different sizes:
        circuit = make_circuit(lambda x,y:  x>y[0:4] , 5)        
        assert( circuit.encrypt_run_decrypt(seq4, seq4) )

        ## < operand
        circuit = make_circuit(lambda x,y:  x<y , 5)

        # close sequences:
        assert( not circuit.encrypt_run_decrypt(seq1, seq3) )        
        # identical sequences
        assert( not circuit.encrypt_run_decrypt(seq1, seq1) )
        # sequences with different sizes:
        circuit = make_circuit(lambda x,y:  x[0:4]<y , 5)        
        assert( circuit.encrypt_run_decrypt(seq4, seq4) )        

        ## >= operand
        circuit = make_circuit(lambda x,y:  x>=y , 5)

        # close sequences:
        assert( circuit.encrypt_run_decrypt(seq1, seq3) )        

        ## len operand
        circuit = make_circuit( lambda x,y: fhe.ones(1) if len( x)==5 else fhe.zeros(1), 5 )
        assert( np.all(circuit.encrypt_run_decrypt(seq1, seq2) == np.array([1])) )

        ## single getitem
        def getitem(x,y):
            out=fhe.zeros(5)
            eseq =  x
            for i in range(5):
                out[i]=eseq[i]
            return out
        circuit = make_circuit(getitem, 5, 4)
        assert( circuit.encrypt_run_decrypt(seq1, seq1) == seq1 )

        ## single encrypted getitem
        def getitem_enc(x,y):
            return  x[y[1]]
        circuit = make_circuit(getitem_enc, 5)
        seq_abcde = Seq('ABCDE')
        fhe_seq_maker = circuit._fhe_seq_maker
        assert( circuit.encrypt_run_decrypt(seq1, seq_abcde) == fhe_seq_maker._letters_to_integers[seq1[fhe_seq_maker._letters_to_integers[seq_abcde[1]]]] )

        ## multiple getitem
        def getitems(x,y):
            eseq =  x
            return eseq[0:2].to_array()
        circuit = make_circuit(getitems, 5, True)
        assert( circuit.encrypt_run_decrypt(seq1, seq1) == seq1[0:2] )

        ## multiple getitem with copy
        def getitemscp(x,y):
            eseq = x
            eseq2= eseq[0:2]
            eseq2[0] = eseq2[1]
            return eseq[0:2].to_array()
        circuit = make_circuit(getitemscp, 5, True)
        assert( circuit.encrypt_run_decrypt(seq1, seq1) == seq1[0:2] )        

        # add sequences (concat them)
        circuit = make_circuit( lambda x,y: (x+y).to_array() , 5, True)
        assert( circuit.encrypt_run_decrypt(seq1, seq2) == (seq1+seq2) )

        ## >= operand for 1 letter
        seq1_2 = Seq('C')
        seq2_2 = Seq('A')

        circuit2 =  make_circuit(lambda x,y:  x>=y , 1)
        assert( circuit2.encrypt_run_decrypt(seq1_2, seq2_2) )
        assert( circuit2.encrypt_run_decrypt(seq1_2, seq1_2) )
        assert( not circuit2.encrypt_run_decrypt(seq2_2, seq1_2) )

    def test_iter(self):
        seq1 = Seq('ACGT*')
        seq2 = Seq('AACGT')
       
        def iter_(x,y):
            seq= x
            seq2=FheMutableSeq(seq[0], x._fhe_seq_maker)
            for c in seq:
                seq2.append(c)
            return seq2.to_array()
        circuit = make_circuit(iter_, 5, True ) 
        assert( circuit.encrypt_run_decrypt(seq1, seq1) == seq2+Seq("*")) 

    def test_startswith(self):
        seq1 = Seq('ACGT')        
        seq2 = Seq('CGTA')
        circuit = make_circuit(lambda x,y: x.startswith(y[0:2]), 4)
        assert( circuit.encrypt_run_decrypt(seq1, seq1) )
        assert( not circuit.encrypt_run_decrypt(seq1, seq2) )

        # also test with naked array:
        circuit = make_circuit(lambda x,y: x.startswith(y[0:2].to_array()), 4)
        assert( circuit.encrypt_run_decrypt(seq1, seq1) )

        # also test with start or end
        circuit = make_circuit(lambda x,y: x.startswith(y[0:2],1), 4)
        assert( circuit.encrypt_run_decrypt(seq1, seq2) )
        circuit = make_circuit(lambda x,y: x.startswith(y[3],None, 3), 4)
        assert( circuit.encrypt_run_decrypt(seq1, seq2) )   

    def test_endswith(self):
        seq1 = Seq('ACGT')        
        seq2 = Seq('CGTA')
        circuit = make_circuit(lambda x,y: x.endswith(y[1:3]), 4)
        assert( not circuit.encrypt_run_decrypt(seq1, seq1) )
        assert( circuit.encrypt_run_decrypt(seq1, seq2) )

    def test_complement(self):
        seq1 = Seq('ACGT')        
        circuit = make_circuit(lambda x,y: x.complement().to_array(), 4, True)
        assert( circuit.encrypt_run_decrypt(seq1, seq1) == seq1.complement())

        # test changing mutable seq in place
        seq2 = MutableSeq('ACGT')
        circuit = make_circuit(lambda x,y:x.complement(True).to_array(), 4, True)
        assert( circuit.encrypt_run_decrypt(seq2, seq2) == seq2.complement(True))

    def test_complement_rna(self):
        seq1 = Seq('ACGU')        
        circuit = make_circuit(lambda x,y: x.complement_rna().to_array(), 4, True)
        assert( circuit.encrypt_run_decrypt(seq1, seq1) == seq1.complement_rna())

        # test changing mutable seq in place
        seq2 = MutableSeq('ACGU')
        circuit = make_circuit(lambda x,y:x.complement(True).to_array(), 4, True)
        assert( circuit.encrypt_run_decrypt(seq2, seq2) == seq2.complement(True))            

    def test_reverse_complement(self):
        seq1 = Seq('ACGT')        
        circuit = make_circuit(lambda x,y: x.reverse_complement().to_array(), 4, True)
        assert( circuit.encrypt_run_decrypt(seq1, seq1) == seq1.reverse_complement())

    def test_reverse_complement_rna(self):
        seq1 = Seq('ACGU')        
        circuit = make_circuit(lambda x,y: x.reverse_complement_rna().to_array(), 4, True)
        assert( circuit.encrypt_run_decrypt(seq1, seq1) == seq1.reverse_complement_rna())

    def test_transcribe(self):
        seq1 = Seq('ACGT')        
        circuit = make_circuit(lambda x,y: x.transcribe().to_array(), 4, True)
        assert( circuit.encrypt_run_decrypt(seq1, seq1) == seq1.transcribe())

    def test_back_transcribe(self):
        seq1 = Seq('ACGU')        
        circuit = make_circuit(lambda x,y: x.back_transcribe().to_array(), 4, True)
        assert( circuit.encrypt_run_decrypt(seq1, seq1) == seq1.back_transcribe())

    def test_translate(self):
        seq1 = Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')        
        circuit = make_circuit(lambda x,y: x.translate().to_array(), len(seq1), True)
        assert( circuit.encrypt_run_decrypt(seq1, seq1) == seq1.translate())

    def test_join(self):
        seq1 = Seq('ACG')
        seq2 = Seq('TTT')
       
        # test joining a FheSeq
        circuit = make_circuit(lambda x,y: x.join(y).to_array(), 3, True)
        assert( circuit.encrypt_run_decrypt(seq2, seq1) == seq2.join(seq1))
        
        # test joining an array
        circuit = make_circuit(lambda x,y: x.join(y.to_array()).to_array(), 3, True)
        assert( circuit.encrypt_run_decrypt(seq2, seq1) == seq2.join(seq1))
        
        # test joining a sequence of FheSeq and arrays
        circuit = make_circuit(lambda x,y: x.join([x,x,y.to_array()]).to_array(), 3, True)
        assert( circuit.encrypt_run_decrypt(seq2, seq1) == seq2.join([seq2,seq2,seq1]))   

    def test_setitem(self):
        seq1 = Seq('XXXXXXXX')
        seq2 = Seq('ABCDEFGH')
       
        def setitem(x,y):
            seq=x
            seq[0]=y[0]
            seq[2:4]=y[2:4]
            seq[5:]=y[5:]
            return seq.to_array()
        circuit = make_circuit( setitem, 8, True )
        assert( circuit.encrypt_run_decrypt(seq1, seq2) == Seq('AXCDXFGH'))

    def test_delitem(self):
        seq1 = Seq('ABCDEFGH')
        seq2 = MutableSeq('ABCDEFGH')
       
        def delitem(x,y):
            seq=x
            del seq[0]
            del seq[2]
            del seq[2:4]
            return seq.to_array()
        circuit = make_circuit( delitem , 8, True)
        del seq2[0]
        del seq2[2]
        del seq2[2:4]
        assert( circuit.encrypt_run_decrypt(seq1, seq1) == seq2)

    def test_append(self):
        seq1 = Seq('AC')
        seq1_bis = MutableSeq('AC')
        seq2 = MutableSeq('GT')
       
        def append(x,y):
            seq=x
            seq2=y
            seq.append(seq2[0])
            return seq.to_array()
        circuit = make_circuit( append, 2, True ) 
        seq1_bis.append(seq2[0])
        assert( circuit.encrypt_run_decrypt(seq1, seq2) == seq1_bis)  

    def test_insert(self):
        seq1 = Seq('AC')
        seq1_bis = MutableSeq('AC')
        seq2 = MutableSeq('GT')
       
        def insert(x,y):
            seq=x
            seq2=y
            seq.insert(1,seq2[0])
            return seq.to_array()
        circuit = make_circuit( insert, 2, True ) 
        seq1_bis.insert(1,seq2[0])
        assert( circuit.encrypt_run_decrypt(seq1, seq2) == seq1_bis)  

    def test_pop(self):
        seq1 = Seq('ACGT')
        seq1_bis = MutableSeq('ACGT')
       
        def pop(x,y):
            seq=x
            seq.pop()
            seq.pop(1)
            return seq.to_array()
        circuit = make_circuit( pop, 4, True ) 
        seq1_bis.pop()
        seq1_bis.pop(1)
        assert( circuit.encrypt_run_decrypt(seq1, seq1) == seq1_bis)  

    def test_extend(self):
        seq1 = Seq('AC')
        seq1_bis = MutableSeq('AC')
        seq2 = MutableSeq('GT')
       
        def extend(x,y):
            seq=x
            seq.extend(y)
            seq.extend(y)
            return seq.to_array()
        circuit = make_circuit( extend, 2, True ) 
        seq1_bis.extend(seq2)
        seq1_bis.extend(seq2)
        assert( circuit.encrypt_run_decrypt(seq1, seq2) == seq1_bis)   



unittest.main()

# suite = unittest.TestLoader().loadTestsFromName('test_FheSeq.TestFheSeq.test_endswith')
# unittest.TextTestRunner(verbosity=1).run(suite)
