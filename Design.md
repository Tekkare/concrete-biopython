# About the Design Choices in concrete.biopython

Let's briefly elucidate the design decisions behind the creation of the **concrete.biopython** library.

## Objective: Alignment with Biopython
The primary goal is to align the `concrete.biopython` library with the principles adopted by the `concrete.fhe` Python library. In the case of `concrete.fhe`, the approach is to ensure compatibility and interoperability with the `numpy` library. This is achieved by offering a compatibility layer between the `concrete.fhe.array` class and `numpy.array` class, allowing seamless transition between unencrypted and encrypted states using `numpy` functions with `concrete.fhe.array` objects.

In the context of `concrete.fhe`, the native objects that undergo encryption and decryption are integer arrays, facilitating a smooth transition from unencrypted to encrypted states. However, the `Biopython` and `concrete.biopython` classes cannot undergo encryption and decryption as straightforwardly.

### Dual Structure with `Seq` / `FheSeq`
The chosen design approach, which has been implemented, revolves around ensuring seamless compatibility with `Biopython` on the unencrypted side, while leveraging a mirrored FHE version exclusively on the encrypted side. The unencrypted side makes use of `Biopython` for tasks like file loading, data preprocessing, and handling `Seq` sequence objects. Subsequently, these objects are transformed into arrays, encrypted, and used to initialize `concrete.biopython` objects such as `FheSeq`, which represent encrypted sequences. After FHE processing, a reverse process is employed for decryption.

The key advantage of this approach lies in the clear distinction between unencrypted and encrypted processing, allowing for the maximum utilization of `Biopython` in the unencrypted part of the code.  
What's more, the BioCircuit circuit wrapper can be used in a majority of usecases (when all inputs are sequences) and facilitates greatly the interfacing between `Seq`
 and `FheSeq`.