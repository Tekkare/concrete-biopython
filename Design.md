# About the Design Choices in concrete.biopython

Let's briefly elucidate the design decisions behind the creation of the **concrete.biopython** library.

## Objective: Alignment with Biopython
The primary goal is to align the `concrete.biopython` library with the principles adopted by the `concrete.fhe` Python library. In the case of `concrete.fhe`, the approach is to ensure compatibility and interoperability with the `numpy` library. This is achieved by offering analogous functions in the `concrete.fhe.array` class as those in the `numpy.array` class, allowing seamless transition between unencrypted and encrypted states using `numpy` functions with `concrete.fhe.array` objects.

In the context of `concrete.fhe`, the native objects that undergo encryption and decryption are integer arrays, facilitating a smooth transition from unencrypted to encrypted states. However, the `Biopython` and `concrete.biopython` classes cannot undergo encryption and decryption as straightforwardly. Consequently, two possible structural approaches emerge.

### I. First Option (Implemented): Dual Structure with `Seq` / `FheSeq`
The chosen design approach, which has been implemented, revolves around ensuring seamless compatibility with `Biopython` on the unencrypted side, while leveraging a mirrored FHE version exclusively on the encrypted side. The unencrypted side makes use of `Biopython` for tasks like file loading, data preprocessing, and handling `Seq` sequence objects. Subsequently, these objects are transformed into arrays, encrypted, and used to initialize `concrete.biopython` objects such as `FheSeq`, which represent encrypted sequences. After FHE processing, a reverse process is employed for decryption.

The key advantage of this approach lies in the clear distinction between unencrypted and encrypted processing, allowing for the maximum utilization of `Biopython` in the unencrypted part of the code.  
However, a drawback is the need for manual object conversions to and from arrays when transitioning between unencrypted and encrypted states, accomplished using the `concrete.biopython.SeqWrapper` class, as detailed in the tutorial.

### II. Second option: `concrete.biopython.Seq` Wrapper Class
An alternative possibility, currently found in the `Seq` branch of the repository, specifically in `tutorial_bis.ypinb`, involves extending the `concrete.biopython` library with an internally concealed wrapper for `Biopython`. This enhancement aims to make `concrete.biopython` objects usable much like `Biopython` in unencrypted processing, while retaining the ability to handle encrypted processing. Under this scenario, a wrapper for `fhe.Compiler` could be established within `concrete.biopython.Seq.Compiler`, allowing for direct encryption and decryption of `Seq` objects, with the translation to and from integer arrays hidden internally.

One drawback of this approach is that the distinction between encrypted and unencrypted objects may become less evident in the code, as well as the differentiation between the capabilities of `concrete.biopython` and `Biopython` objects.  
Nevertheless, the benefit lies in a more direct and straightforward method for encrypting and decrypting `Seq` objects (and potentially other class objects in the future).

### III. Third option: hybrid solution

A third option could be to keep the separation suggested in the first solution, while using a compiler wrapper allowing to easily compile functions onto **Seq** sequences objects.


## Features

- FheSeq functions: the implemented (and todo) functions of **FheSeq** classes are the same as those of **Biopython.Seq** classes (when they are possible to do in FHE).
- Alphabet setting for bitwidth optimization: the FheSeq object use a default alphabet containing all letters, but it can be set to fewer specific values which improves the FHE speed.