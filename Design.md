# About the concrete.biopython design choices
Let's explain breifly the design choices made to build the **concrete.biopython** library.

## Objective: aligning the lib with Biopython
The objective is to follow what has been done with the concrete python **concrete.fhe** library. This library aligns with the library **numpy** so as to be maximally compatible (same functions in **concrete.fhe.array** class than in **numpy.array** class) and interoperable (using **numpy** functions on **concrete.fhe.array** objects).  

In **concrete.fhe** the native objects that are encrypted and decrypted are integer arrays, which allows a smooth usage of **concrete.fhe.array** objects from unencrypted to encrypted states. However, the **Biopython** and **concrete.biopython** classes cannot be encrypted and decrypted in a straightforward way. That's why two possible structure are possible.

## First option (implemented): Dual structure with Seq / FheSeq
This is the design choice which was implementend. The idea is to keep everything working under **Biopython** on the unencrypted side, and to use a mirror FHE version only on the encrypted side. The unencrypted side thus uses **Biopython** for loading files, preprocesssing data, and loading **Seq** sequence objects. Then, objects are transformed to arrays, encrypted, and the encrypted arrays are used to initialize **concrete.biopython** objects such as **FheSeq** encrypted sequences. After FHE processing, the reverse process is used for decryption.  

The main advantage here is that there is a clear distinction between unencrypted and encrypted processing, with maximum and direct use of **Biopython** on the unencrypted side.  
The inconvenient is that the switching from unencrypted to encrypted (and vice-versa) requires manual conversions of objects to arrays and back to objects, using the **concrete.biopython.SeqWrapper** class (as explained in the tutorial).

## Second option : concrete.biopython.Seq wrapper class

A second possibility that is currently presented in the repo's branch **Seq**, in `tutorial_bis.ypinb`, is to extend the library **concrete.biopython** with a hidden inner wrapping of **Biopython**, so as to make the objects of **concrete.biopython** usable like **Biopython** in unencrypted processing, in addition to being able to deal with encrypted processing. In this case, we could create a wrapper for **fhe.Compiler** in **concrete.biopython.Seq.Compiler** that would directly encrypt and decrypt **Seq** objects, and hide the in and out translation to integer arrays inside.

The inconvenient is that the distinction between encrypted and unencrypted objects would be less clear in the code, as well as the distinction between what can **concrete.biopython** and **Biopython** object do or not.
The advantage is a more direct and straightfoward way to encrypt and decrypt **Seq** objects (and later other class objects).



# About the Design Choices in concrete.biopython

Let's briefly elucidate the design decisions behind the creation of the **concrete.biopython** library.

## Objective: Alignment with Biopython
The primary goal is to align the `concrete.biopython` library with the principles adopted by the `concrete.fhe` Python library. In the case of `concrete.fhe`, the approach is to ensure compatibility and interoperability with the `numpy` library. This is achieved by offering analogous functions in the `concrete.fhe.array` class as those in the `numpy.array` class, allowing seamless transition between unencrypted and encrypted states using `numpy` functions with `concrete.fhe.array` objects.

In the context of `concrete.fhe`, the native objects that undergo encryption and decryption are integer arrays, facilitating a smooth transition from unencrypted to encrypted states. However, the `Biopython` and `concrete.biopython` classes cannot undergo encryption and decryption as straightforwardly. Consequently, two possible structural approaches emerge.

## First Option (Implemented): Dual Structure with `Seq` / `FheSeq`
The chosen design approach, which has been implemented, revolves around ensuring seamless compatibility with `Biopython` on the unencrypted side, while leveraging a mirrored FHE version exclusively on the encrypted side. The unencrypted side makes use of `Biopython` for tasks like file loading, data preprocessing, and handling `Seq` sequence objects. Subsequently, these objects are transformed into arrays, encrypted, and used to initialize `concrete.biopython` objects such as `FheSeq`, which represent encrypted sequences. After FHE processing, a reverse process is employed for decryption.

The key advantage of this approach lies in the clear distinction between unencrypted and encrypted processing, allowing for the maximum utilization of `Biopython` in the unencrypted part of the code.  
However, a drawback is the need for manual object conversions to and from arrays when transitioning between unencrypted and encrypted states, accomplished using the `concrete.biopython.SeqWrapper` class, as detailed in the tutorial.

## Second option: `concrete.biopython.Seq` Wrapper Class
An alternative possibility, currently found in the `Seq` branch of the repository, specifically in `tutorial_bis.ypinb`, involves extending the `concrete.biopython` library with an internally concealed wrapper for `Biopython`. This enhancement aims to make `concrete.biopython` objects usable much like `Biopython` in unencrypted processing, while retaining the ability to handle encrypted processing. Under this scenario, a wrapper for `fhe.Compiler` could be established within `concrete.biopython.Seq.Compiler`, allowing for direct encryption and decryption of `Seq` objects, with the translation to and from integer arrays hidden internally.

One drawback of this approach is that the distinction between encrypted and unencrypted objects may become less evident in the code, as well as the differentiation between the capabilities of `concrete.biopython` and `Biopython` objects.  
Nevertheless, the benefit lies in a more direct and straightforward method for encrypting and decrypting `Seq` objects (and potentially other class objects in the future).

## Third option: hybrid solution

A third option could be to keep the separation suggested in the first solution, while using a compiler wrapper allowing to easily compile functions onto **Seq** sequences objects.
