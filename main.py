from concrete import fhe
import numpy as np

from BioConcrete.FheSeq import FheSeq
from Bio.Seq import Seq, MutableSeq
from BioConcrete.utils import seqToIntegers, seqFromIntegers



def dna_complement(encrypted_integers):
    encSeq = FheSeq(encrypted_integers)
    return encSeq.complement().toArray()

def dna_reverse_complement(encrypted_integers):
    encSeq = FheSeq(encrypted_integers)
    return encSeq.reverse_complement().toArray()


class BioConcreteCircuit:
   
    def __init__(self, seq_length):
        self.seq_length=seq_length
        self.inputset=[
            np.random.randint(0, 4, size=(seq_length,))
            for _ in range(100)
        ]
        self.circuit = None

    def set(self, circuitFunction, verbose=False):
        compiler = fhe.Compiler(lambda data: circuitFunction(data), {"data": "encrypted"})
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
    
    def run(self, seq, simulate=False):
        if not self.circuit: raise Error('circuit was not set')
        assert len(seq) == self.seq_length, f"Sequence length is not correct, should be {self.seq_length} characters"

        # convert letters to integers
        integers = seqToIntegers(seq)
        output_seq = self.circuit.simulate(integers) if simulate else self.circuit.encrypt_run_decrypt(integers)

        # convert back integers to letters    
        return seqFromIntegers(output_seq)



import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--seq', type=str, help='DNA sequence')
parser.add_argument('--sim', action='store_true', help='simulate encryption')
args = parser.parse_args()

if not args.seq:
    args.seq = 'ACGT'


circuit = BioConcreteCircuit(len(args.seq))
seq = Seq(args.seq)
print('processing sequence : ', seq)


#### complement
print('\ncomplement :')

circuit.set(dna_complement)
decrypted_seq = circuit.run(seq, args.sim)

print("decrypted :", decrypted_seq)
print("expected :", seq.complement())


#### reverse complement
print('\nreverse complement :')

circuit.set(dna_reverse_complement)
decrypted_seq = circuit.run(seq, args.sim)
print("decrypted :", decrypted_seq)
print("expected :", seq.reverse_complement())
