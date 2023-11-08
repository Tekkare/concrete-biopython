import numpy as np
from concrete import fhe

from Bio.Seq import Seq, MutableSeq

from concrete_biopython.FheSeq import FheSeq, FheMutableSeq
from concrete_biopython.SeqWrapper import SeqWrapper


# # Separate sequences in two groups 
# def separate_sequences(seq_list, encryption):
#     # Check if the lengths of the list and dictionary match
#     if len(seq_list) != len(encryption):
#         raise ValueError("Sequence list and encryption dictionary sizes do not match")

#     # Initialize two lists to store clear and encrypted values
#     clear_values = []
#     encrypted_values = []

#     # Iterate over the elements of the list and the dictionary
#     for i in range(len(seq_list)):
#         value = encryption[key]
#         if value == "clear":
#             clear_values.append(seq_list[i])
#         elif value == "encrypted":
#             encrypted_values.append(seq_list[i])
#         else:
#             raise ValueError(f"Invalid value in the dictionary for key {key}: {value}\n\
#                 should be either \"encrypted\" or \"clear\"")

#     return encrypted_values, clear_values

# def concat_sequences(seq_list):
#     """
#     Concatenate a list of Seq or FheSeq sequences, and records their indices as slice objects.
#     """
#     seq= seq_list[0]
#     index=len(seq_list[0])
#     slices=[slice(0,index)]
#     for i in range(1, len(seq_list)):
#         seq += seq_list[i] # Seq objects can be directly added
#         slices.append(slice(index, index+len(seq_list[i])))
#         index+=len(seq_list[i])
#     return seq, slices

# def slice_sequence(seq, slices):
#     """
#     Slice a long sequence into its original sub-sequences
#     """
#     seq_list=[]
#     index=0
#     for ind in slices:
#         seq_list.append(seq[ind.start:ind.stop])
#     return seq_list

# # Compute the minimal alphabet to represent a Seq object
# def compute_alphabet(seq):
#     # Use a set to store unique letters
#     unique_letters = set(str(seq))
#     # Convert the set of unique letters back to a sorted string
#     return ''.join(sorted(unique_letters))


# Compute the minimal alphabet to represent a Seq object
def compute_alphabet(seq_list):
    # Use a set to store unique letters
    unique_letters = set()
    for seq in seq_list:
        # update with every seq letters
        unique_letters.update(str(seq))
    # Convert the set of unique letters back to a sorted string
    return ''.join(sorted(unique_letters))


# Wrap a function so that it can use MutableSeq objects in inputs and a MutableSeq/Seq in output
def function_wrapper_factory(function, param_names):

    # Create a dynamic version of the above function wrapper
    # than can take named parameters, the ones from encryption:
    #
    # def function_wrapper(seq1, seq2, ...):
    #     fhe_arrays = [seq1, seq2, ...]
    #     fhe_mutableSeqs = [ FheMutableSeq(fhe_seq) for fhe_seq in fhe_arrays ]
    #     output = function(*fhe_mutableSeqs)
    #     # if the output of the function is a FheSeq or FheMutableSeq, convert it back to an array
    #     if isinstance(output, FheMutableSeq) or isinstance(output, FheSeq):
    #         output = output.toArray()
    #     return output    

    # Create the function definition as a string
    func_str = f"def function_wrapper({', '.join(param_names)}):\n"+
        f"    fhe_arrays = [{', '.join(param_names)}]\n"+
        "    fhe_mutableSeqs = [ FheMutableSeq(fhe_seq) for fhe_seq in fhe_arrays ]\n"+
        "    output = function(*fhe_mutableSeqs)\n"+
        "    if isinstance(output, FheMutableSeq) or isinstance(output, FheSeq):\n"+
        "       output = output.toArray()\n"+
        "    return output"

    # Compile the function string
    compiled_func = compile(func_str, '<string>', 'exec')
    
    # Create a namespace for the function
    local_vars = {}
    
    # Execute the compiled function in the local namespace
    exec(compiled_func, globals(), local_vars)
    
    # Get a reference to the dynamically created function
    dynamic_function_wrapper = local_vars['function_wrapper']
    
    return dynamic_function_wrapper
    

class BioCircuit:
    """
    A wrapper class of concrete circuit that uses Seq objects as inputs and outputs
    """

    def __init__(self,
                 function,
                 encryption,
                 seq_list,
                 configuration,
                 **kwargs
                ):
        """
        Arguments for the constructor:
        
        - `function` a function taking FheMutableSeq objects
        - `encryption` a dictionnary with string keys and values equal to
          "encrypted" or "clear" as required to compile a circuit
        - `seq_list` a list of Seq inputs that can be used by the circuit, this is
          required to compute the alphabet and the Seq input lenghts
        - `configuration` the configuration for the compiler
        - `**kwargs` any other named arguments for the compiler
        """

        # verify the seq_list types
        for seq in seq_list:
            if not isinstance(seq, Seq) or isinstance(seq, MutableSeq):
                raise ValueError("Sequence must be of type Seq or MutableSeq")        

        # verify the encryption
        for key in encryption.keys():
            value = encryption[key]
            if value != "clear" and value != "encrypted":
                raise ValueError(f"Invalid value in the dictionary for key {key}: {value}\n\
                    should be either \"encrypted\" or \"clear\"")    

        # # concatenate the seq_list into one long Seq
        # seq, _ = concat_sequences(seq_list)
        # # compute the minimal alphabet to represent the seq_list
        # self._alphabet = compute_alphabet(seq)

        # compute the minimal alphabet to represent the seq_list
        self._alphabet = compute_alphabet(seq_list)
        # set the alphabet to be able to call SeqWrapper.maxInteger()
        SeqWrapper.setAlphabet(self._alphabet)

        # Create the circuit integer seq_list from the seq
        # Use SeqWrapper.maxInteger() to know the maximum integer that can be
        # used to represent a character in FheSeq obects        
        inputset = [np.random.randint(0, SeqWrapper.maxInteger()+1, size=(len(seq),))
                    for _ in range (300)]

        # Wrap the function so that it can use any number of MutableSeq objects in inputs and a MutableSeq/Seq in output
        # get parameter names from encryption directory
        param_names = list(encryption.keys())
        wrapped_function = function_wrapper_factory(function, param_names)

        # compile the circuit   
        compiler = fhe.Compiler(wrapped_function, encryption)
        self._circuit = compiler.compile(
            inputset = inputset,
            configuration = configuration,
            **kwargs,
        )
        # reset the alphabet in case we use another circuit before running this one
        SeqWrapper.resetAlphabet()

    def encrypt(self, *seq_list):
        # convert sequences to integer representation
        integer_list = [];
        for seq in seq_list:
            if not isinstance(seq, Seq) or isinstance(seq, MutableSeq):
                raise ValueError("All sequences must be of type Bio.Seq.Seq or Bio.Seq.MutableSeq")
            integer_list.append( SeqWrapper.toIntegers(seq) )
        # return the list of sequences as integers encrypted in FHE
        encrypted_input = self._circuit.encrypt(*integer_list)
        return encrypted_input

    def run(self, encrypted_input):
        # set the alphabet for improved computation performances
        SeqWrapper.setAlphabet(self._alphabet)
        encrypted_output = self._circuit.run(encrypted_input)
        # always reset the alphabet after running the circuit
        SeqWrapper.resetAlphabet()
        return encrypted_output

    def decrypt(self, encrypted_output):
        return self._circuit.decrypt(encrypted_output)

    def encrypt_run_decrypt(self, *seq_list):
        encrypted_input = self._circuit.encrypt(*seq_list)
        encrypted_output = self._circuit.run(encrypted_input)
        decrypted_result = self._circuit.decrypt(encrypted_output)
        return decrypted_result


