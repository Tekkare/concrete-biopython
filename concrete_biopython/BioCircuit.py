import numpy as np
from concrete import fhe

from Bio.Seq import Seq, MutableSeq

from concrete_biopython.FheSeq import FheSeq, FheMutableSeq
from concrete_biopython.SeqWrapper import SeqWrapper


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
def function_wrapper_factory(param_names):

    """
    Create a dynamic version of the above double function wrapper
    than can take named parameters, the ones from encryption:
    
    def function_double_wrapper(function):
    
        def function_wrapper(seq1, seq2, ...):
            fhe_arrays = [seq1, seq2, ...]
            fhe_mutableSeqs = [ FheMutableSeq(fhe_seq) for fhe_seq in fhe_arrays ]
            output = function(*fhe_mutableSeqs)
            # if the output of the function is a FheSeq or FheMutableSeq, convert it back to an array
            if isinstance(output, FheMutableSeq) or isinstance(output, FheSeq):
                output = output.toArray()
            return output  
    
        return function_wrapper  
    """

    # Create the function definition as a string
    func_str =  "def function_double_wrapper(function):\n"
    func_str += f"    def function_wrapper({', '.join(param_names)}):\n"
    func_str += f"        fhe_arrays = [{', '.join(param_names)}]\n"
    func_str +=  "        fhe_mutableSeqs = [ FheMutableSeq(fhe_seq) for fhe_seq in fhe_arrays ]\n"
    func_str +=  "        output = function(*fhe_mutableSeqs)\n"
    func_str +=  "        if isinstance(output, FheMutableSeq) or isinstance(output, FheSeq):\n"
    func_str +=  "            output = output.toArray()\n"
    func_str +=  "        return output\n"
    func_str +=  "    return function_wrapper"

    # Compile the function string
    compiled_func = compile(func_str, '<string>', 'exec')
    
    # Create a namespace for the function
    local_vars = {}
    
    # Execute the compiled function in the local namespace
    exec(compiled_func, globals(), local_vars)
    
    # Get a reference to the dynamically created function
    dynamic_function_wrapper = local_vars['function_double_wrapper']
    
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
        #SeqWrapper.setAlphabet(self._alphabet)

        # Create the circuit integer seq_list from the seq
        # Use SeqWrapper.maxInteger() to know the maximum integer that can be
        # used to represent a character in FheSeq obects        
        inputset = [tuple([np.random.randint(0, SeqWrapper.maxInteger()+1, size=(len(seq),))
                            for seq in seq_list]) for _ in range (300)]

        # Wrap the function so that it can use any number of MutableSeq objects in inputs and a MutableSeq/Seq in output
        # get parameter names from encryption directory
        param_names = list(encryption.keys())
        wrapped_function = function_wrapper_factory(param_names)(function)

        # compile the circuit   
        compiler = fhe.Compiler(wrapped_function, encryption)    
        self._circuit = compiler.compile(
            inputset = inputset,
            configuration = configuration,
            **kwargs,
        )        
        # reset the alphabet in case we use another circuit before running this one
        #SeqWrapper.resetAlphabet()

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


