import numpy as np
from concrete import fhe

import time

from Bio.Seq import Seq, MutableSeq

from concrete_biopython.FheSeq import SeqInterface, FheSeq, FheMutableSeq


def measure_time(function, descripton, *inputs):
    """
    Compute a function on inputs and return output along with duration
    """
    start = time.time()
    output = function(*inputs)
    end = time.time()
    print(f"|  {descripton} : {end-start:.2f} s  |")
    return output


def function_double_wrapper_factory(param_names):
    """
    Our goal is to get a function that can use MutableSeq objects in inputs and a MutableSeq/Seq in output
    As the number of input arrays is variable, and as the compiler needs a function with names to compile,
    we need to create this function dynamically during the execution, when we know the names and number of
    array parameters.
    To this purpose, we will dynamically create a function function_double_wrapper from string, that given
    a function process_sequence (the sequence processing function) as parameter, will output a function function_wrapper
    taking the correct number of input arrays, convert them to MutableSeq objects, and apply process_sequence to them
    """

    """
    Create a dynamic version of the above double function wrapper that can take named parameters (the ones from encryption),
    and output a function wrapper that takes a function process_sequence, which will in turn apply a function process_sequence
    onto FheMutableSeq objects made from the arrays in parameter
    
    def function_double_wrapper(process_sequence, seq_interface):
    
        def function_wrapper(seq1, seq2, ...):
            fhe_arrays = [seq1, seq2, ...]
            fhe_mutableSeqs = [ seq_interface.FheMutableSeq(fhe_arr) for fhe_arr in fhe_arrays ]
            output = process_sequence(*fhe_mutableSeqs)
            # if the output of the function is a FheSeq or FheMutableSeq, convert it back to an array
            if isinstance(output, FheMutableSeq) or isinstance(output, FheSeq):
                output = output.to_array()
            if isinstance(output, list):
                for el in output:
                    el = el.to_array() if isinstance(el, FheMutableSeq) or isinstance(el, FheSeq) else el
            return output  
    
        return function_wrapper
    """

    # Create the function definition as a string
    func_str = "def function_double_wrapper(process_sequence, seq_interface):\n"
    func_str += f"    def function_wrapper({', '.join(param_names)}):\n"
    func_str += f"        fhe_arrays = [{', '.join(param_names)}]\n"
    func_str += "        fhe_mutableSeqs = [ seq_interface.FheMutableSeq(fhe_arr) for fhe_arr in fhe_arrays ]\n"
    func_str += "        output = process_sequence(*fhe_mutableSeqs)\n"
    func_str += (
        "        if isinstance(output, FheMutableSeq) or isinstance(output, FheSeq):\n"
    )
    func_str += "            output = output.to_array()\n"
    func_str += "        if isinstance(output, list) or isinstance(output, tuple):\n"
    func_str += "            for i in range(0,len(output)):\n"
    func_str += "                el = output[i]\n"
    func_str += "                output[i] = el.to_array() if isinstance(el, FheMutableSeq) or isinstance(el, FheSeq) else el\n"
    func_str += "        return output\n"
    func_str += "    return function_wrapper"

    # Compile the function string
    compiled_func = compile(func_str, "<string>", "exec")

    # Create a namespace for the function
    local_vars = {}

    # Execute the compiled function in the local namespace
    exec(compiled_func, globals(), local_vars)

    # Get a reference to the dynamically created function
    dynamic_function_wrapper = local_vars["function_double_wrapper"]

    return dynamic_function_wrapper


class BioCircuit:
    """
    A wrapper class of concrete circuit that uses Seq objects as inputs and outputs

    TODO:
    - add production APIs (e.g., saving and loading BioCircuits, key management, deployment pipeline)
    - add the possibility to provide an encryption dictionnary telling whether Sequences are encrypted or clear
    for instance: {"seq1": "encrypted", "seq2": "clear"}. This requires concrete to be able to deals with clear Tracers
    """

    def __init__(
        self,
        function,
        # encryption,
        len_seqs,
        seq_interface,
        configuration,
        inputset=None,
        seq_output=False,
        show_timing=True,
        **kwargs,
    ):
        """
        Arguments for the constructor:

        - `function` a function taking FheMutableSeq objects
        - `len_seqs` a list of lengths of the Seq inputs that can be used by the circuit
        - `seq_interface` a SeqInterface object to use for creating FheSeq objects
        - `configuration` the configuration for the compiler
        - `inputset` an optional custom inputset
        - `seq_output` whether outputs are to be converted to Seq objects
        - `show_timing` whether to show timing or not
        - `**kwargs` any other named arguments for the compiler

        TODO : - `encryption` a dictionnary with string keys and values equal to
          "encrypted" or "clear" as required to compile a circuit
        """

        self._len_seqs = len_seqs
        self._show_timing = show_timing
        self._seq_output = seq_output
        self._seq_interface = seq_interface

        assert(isinstance(seq_interface, SeqInterface), "seq_interface must be of type SeqInterface")

        # create an encryption dictionnary from the number of input sequences, where everything is encrypted
        encryption = {"seq" + str(i): "encrypted" for i in range(1, len(len_seqs) + 1)}

        # Create the circuit inputset from the seq lengths and the maximum integer
        # Use seq_interface.max_integer() to know the maximum integer that can be
        # used to represent a character in FheSeq obects
        if inputset is not None:
            # check correctness of inputset
            if not (
                isinstance(inputset, list)
                and all(
                    isinstance(item, tuple)
                    and all(
                        isinstance(seq_item, (Seq, MutableSeq)) for seq_item in item
                    )
                    for item in inputset
                )
            ):
                raise ValueError(
                    "The inputset should be either None or with type List[Tuple[Union[Seq, MutableSeq]]]"
                )

            inputset = [
                tuple([seq_interface.to_integers(seq_item) for seq_item in item])
                for item in inputset
            ]

        else:
            inputset = [
                tuple(
                    [
                        np.random.randint(
                            0, seq_interface.max_integer() + 1, size=(len_seq,)
                        )
                        for len_seq in len_seqs
                    ]
                )
                for _ in range(300)
            ]

        # Wrap the function so that it can use any number of MutableSeq objects in inputs and a MutableSeq/Seq in output
        # get parameter names from encryption directory
        param_names = list(encryption.keys())
        wrapped_function = function_double_wrapper_factory(param_names)(
            function, seq_interface
        )

        # compile the circuit
        compiler = fhe.Compiler(wrapped_function, encryption)

        compile_ = lambda _: compiler.compile(
            inputset=inputset,
            configuration=configuration,
            **kwargs,
        )

        if show_timing:
            self._circuit = measure_time(compile_, "Compiling ", None)
        else:
            self._circuit = compile_(None)

    def encrypt(self, *seq_list):
        """
        Encrypt a list of Seq or MutableSeq objects
        """

        # convert sequences to integer representation
        integer_list = []
        for i in range(0, len(seq_list)):
            seq = seq_list[i]
            if not (isinstance(seq, (Seq, MutableSeq))):
                raise ValueError("All sequences must be of type Seq or MutableSeq")
            if not len(seq) == self._len_seqs[i]:
                raise ValueError(
                    f"Sequence number {i} has not the right length. Expected {self._len_seqs[i]}, got {len(seq)}"
                )
            integer_list.append(self._seq_interface.to_integers(seq))

        if self._show_timing:
            encrypted_input = measure_time(
                self._circuit.encrypt, "Encrypting ", *integer_list
            )
        else:
            encrypted_input = self._circuit.encrypt(*integer_list)

        return encrypted_input

    def run(self, encrypted_input):
        """
        Run the circuit on encrypted input
        """
        if self._show_timing:
            encrypted_output = measure_time(
                self._circuit.run, "Running ", encrypted_input
            )
        else:
            encrypted_output = self._circuit.run(encrypted_input)

        return encrypted_output

    def decrypt(self, encrypted_output, as_seq=False):
        """
        Decrypt an convert back to a sequence if required
        """
        if self._show_timing:
            res = measure_time(self._circuit.decrypt, "Decrypting ", encrypted_output)
        else:
            res = self._circuit.decrypt(encrypted_output)

        # decrypt and convert back to a sequence if required
        if self._seq_output:
            res = self._seq_interface.array_to_seq(res)

        return res

    def encrypt_run_decrypt(self, *seq_list):
        """
        Encrypt, run and decrypt at once
        """
        encrypted_input = self.encrypt(*seq_list)
        encrypted_output = self.run(encrypted_input)
        decrypted_result = self.decrypt(encrypted_output)
        return decrypted_result

    def simulate(self, *seq_list):
        """
        Simulate a run of the circuit
        """
        if self._show_timing:
            encrypted_output = measure_time(
                self._circuit.simulate, "Running ", *seq_list
            )
        else:
            encrypted_output = self._circuit.simulate(*seq_list)

        return encrypted_output
