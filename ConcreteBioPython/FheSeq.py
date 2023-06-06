from abc import ABC
from concrete import fhe
import numpy as np
import numbers
from .SeqWrapper import SeqWrapper


class _SeqAbstractBaseClass(ABC):

    _DNAcomplementTable = fhe.LookupTable(SeqWrapper.get_DNA_complementTable()) 
    _RNAcomplementTable = fhe.LookupTable(SeqWrapper.get_RNA_complementTable()) 
    _transcriptionTable = fhe.LookupTable(SeqWrapper.get_transcriptionTable()) 
    _backTranscriptionTable = fhe.LookupTable(SeqWrapper.get_back_transcriptionTable())
    _translationReductionTable = fhe.LookupTable(SeqWrapper.get_translationReductionTable())

    def __init__(self, data, length=None):
        if data is None:
            if length is None:
                raise ValueError("length must not be None if data is None")
            elif length == 0:
                self._data = fhe.zeros(length)
            elif length < 0:
                raise ValueError("length must not be negative.")
            else:
                self._data = _UndefinedSequenceData(length)
        elif isinstance(data, fhe.tracing.tracer.Tracer):
            self._data = data
        else:
            print(type(data))
            raise TypeError(
                "data should be a concrete.fhe.tracing.tracer.Tracer object"
            )     

    def toArray(self):
        return self._data

    def __eq__(self, other):
        """Compare the sequence to another sequence or an array.

        Sequences are equal to each other if their sequence contents is identical
        """
        if isinstance(other, _SeqAbstractBaseClass):
            return (len(self) - np.sum(self._data == other._data))==0
        elif isinstance(other, fhe.tracing.tracer.Tracer):
            return (len(self) - np.sum(self._data == other))==0
        else:
            return NotImplemented

    # def __lt__(self, other):
    #     """Implement the less-than operand."""
    #     if isinstance(other, _SeqAbstractBaseClass):
    #         return self._data < other._data
    #     elif isinstance(other, fhe.tracing.tracer.Tracer):
    #         return self._data < other
    #     else:
    #         return NotImplemented

    # def __le__(self, other):
    #     """Implement the less-than or equal operand."""
    #     if isinstance(other, _SeqAbstractBaseClass):
    #         return self._data <= other._data
    #     elif isinstance(other, fhe.tracing.tracer.Tracer):
    #         return self._data <= other
    #     else:
    #         return NotImplemented

    # def __gt__(self, other):
    #     """Implement the greater-than operand."""
    #     if isinstance(other, _SeqAbstractBaseClass):
    #         return self._data > other._data
    #     elif isinstance(other, fhe.tracing.tracer.Tracer):
    #         return self._data > other
    #     else:
    #         return NotImplemented

    # def __ge__(self, other):
    #     """Implement the greater-than or equal operand."""
    #     if isinstance(other, _SeqAbstractBaseClass):
    #         return self._data >= other._data
    #     elif isinstance(other, fhe.tracing.tracer.Tracer):
    #         return self._data >= other  
    #     else:
    #         return NotImplemented

    def __len__(self):
        """Return the length of the sequence."""
        return self._data.size

    # tracer.__iter__ is not availabe yet
    # def __iter__(self):
    #     """Return an iterable of the sequence."""
    #     return self._data.__iter__()

    def __getitem__(self, index):
        """Return a subsequence as a single integer or as a sequence object.
        """
        if isinstance(index, numbers.Integral):
            # Return a single integer
            return self._data[index]
        else:
            # Return the (sub)sequence as another Seq/MutableSeq object
            return self.__class__(self._data[index])

    def __add__(self, other):
        """Add a sequence or array to this sequence.
        """
        if isinstance(other, _SeqAbstractBaseClass):
            return self.__class__(np.concatenate((self._data, other._data), axis=0))
        elif isinstance(other, fhe.tracing.tracer.Tracer):
            return self.__class__(np.concatenate((self._data, other), axis=0))
        else:
            return NotImplemented

    # Cannot be implemented, because it collides with array __add__ function            
    # def __radd__(self, other):
    #     """Add an array on the left.
    #     """
    #     if isinstance(other, fhe.tracing.tracer.Tracer):
    #         return self.__class__(fhe.concatenate((other , self._data), axis=0))
    #     else:
    #         return NotImplemented       

    # def __mul__(self, other):
    #     """Multiply sequence by integer.
    #     """
    #     if not isinstance(other, numbers.Integral):
    #         raise TypeError(f"can't multiply {self.__class__.__name__} by non-int type")
    #     # we would like to simply write
    #     # data = self._data * other
    #     # here, but currently that causes a bug on PyPy if self._data is a
    #     # bytearray and other is a numpy integer. Using this workaround:
    #     data = self._data.__mul__(other)
    #     return self.__class__(data)

    # def __rmul__(self, other):
    #     """Multiply integer by sequence.
    #     """
    #     if not isinstance(other, numbers.Integral):
    #         raise TypeError(f"can't multiply {self.__class__.__name__} by non-int type")
    #     # we would like to simply write
    #     # data = self._data * other
    #     # here, but currently that causes a bug on PyPy if self._data is a
    #     # bytearray and other is a numpy integer. Using this workaround:
    #     data = self._data.__mul__(other)
    #     return self.__class__(data)

    # def __imul__(self, other):
    #     """Multiply the sequence object by other and assign.
    #     """
    #     if not isinstance(other, numbers.Integral):
    #         raise TypeError(f"can't multiply {self.__class__.__name__} by non-int type")
    #     # we would like to simply write
    #     # data = self._data * other
    #     # here, but currently that causes a bug on PyPy if self._data is a
    #     # bytearray and other is a numpy integer. Using this workaround:
    #     data = self._data.__mul__(other)
    #     return self.__class__(data)    

    # def count(self, sub, start=None, end=None):
    #     """Return a non-overlapping count, like that of a python string.

    #     The number of occurrences of substring argument sub in the
    #     (sub)sequence given by [start:end] is returned as an integer.
    #     Optional arguments start and end are interpreted as in slice
    #     notation.

    #     Arguments:
    #      - sub - a string or another FheSeq object to look for
    #      - start - optional integer, slice start
    #      - end - optional integer, slice end
    #     """
    #     if isinstance(sub, _SeqAbstractBaseClass):
    #         sub = sub._data
    #     elif isinstance(sub, fhe.tracing.tracer.Tracer):
    #         pass
    #     else:
    #         return NotImplemented

    #     if start:
    #         if end:
    #             sub = sub[start:end]
    #         else:
    #             sub = sub[start:]

    #     count =

    #     return count   

    def complement(self, inplace=False):
        """Return the complement as an RNA sequence.

        Any T in the sequence is treated as a U
        Any other letter than ACGT will be changed to the character with index 0 which should not be done

        The sequence is modified in-place and returned if inplace is True

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if
        ``complement_rna`` is called on a ``Seq`` object with ``inplace=True``.
        """
        complement_seq = fhe.zeros(len(self))
        for i in range(len(self)):
            complement_seq[i] = _SeqAbstractBaseClass._DNAcomplementTable[self._data[i]]
        if inplace:
            if self.__class__ == FheSeq:
                raise TypeError("Sequence is immutable")
            self._data = complement_seq
            return self
        else:
            return self.__class__(complement_seq)

    def complement_rna(self, inplace=False):
        """Return the complement as an RNA sequence.

        Any T in the sequence is treated as a U

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if
        ``complement_rna`` is called on a ``Seq`` object with ``inplace=True``.
        """
        complement_seq = fhe.zeros(len(self))
        for i in range(len(self)):
            complement_seq[i] = _SeqAbstractBaseClass._RNAcomplementTable[self._data[i]]
        if inplace:
            if self.__class__ == FheSeq:
                raise TypeError("Sequence is immutable")
            self._data = complement_seq
            return self
        else:
            return self.__class__(complement_seq)

    def reverse_complement(self, inplace=False):
        """Return the reverse complement as a DNA sequence.

        Any U in the sequence is treated as a T

        The sequence is modified in-place and returned if inplace is True

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if
        ``reverse_complement`` is called on a ``Seq`` object with
        ``inplace=True``.
        """        
        reverse_complement = self.complement(inplace)
        if inplace:
            if not isinstance(self, FheSeq):
                raise TypeError("Sequence is immutable")
        reverse_complement._data[::-1] = reverse_complement._data
        return reverse_complement

    def reverse_complement_rna(self, inplace=False):
        """Return the reverse complement as an RNA sequence.

        The sequence is modified in-place and returned if inplace is True

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if
        ``reverse_complement_rna`` is called on a ``Seq`` object with
        ``inplace=True``.
        """
        reverse_complement_rna = self.complement_rna(inplace)
        if inplace:
            if not isinstance(self, FheSeq):
                raise TypeError("Sequence is immutable")
        reverse_complement_rna._data[::-1] = reverse_complement_rna._data
        return reverse_complement_rna

    def transcribe(self, inplace=False):
        """Transcribe a DNA sequence into RNA and return the RNA sequence as a new Seq object.

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if
        ``transcribe`` is called on a ``Seq`` object with ``inplace=True``.

        Trying to transcribe an RNA sequence has no effect.
        If you have a nucleotide sequence which might be DNA or RNA
        (or even a mixture), calling the transcribe method will ensure
        any T becomes U.

        Trying to transcribe a protein sequence will replace any
        T for Threonine with U for Selenocysteine, which has no
        biologically plausible rational.
        """
        transcribed_seq = fhe.zeros(len(self))
        for i in range(len(self)):
            transcribed_seq[i] = _SeqAbstractBaseClass._transcriptionTable[self._data[i]]
        if inplace:
            if self.__class__ == FheSeq:
                raise TypeError("Sequence is immutable")
            self._data = transcribed_seq
            return self
        else:
            return self.__class__(transcribed_seq)

    def back_transcribe(self, inplace=False):
        """Return the DNA sequence from an RNA sequence by creating a new Seq object.

        The sequence is modified in-place and returned if inplace is True:

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if
        ``transcribe`` is called on a ``Seq`` object with ``inplace=True``.

        Trying to back-transcribe DNA has no effect, If you have a nucleotide
        sequence which might be DNA or RNA (or even a mixture), calling the
        back-transcribe method will ensure any U becomes T.

        Trying to back-transcribe a protein sequence will replace any U for
        Selenocysteine with T for Threonine, which is biologically meaningless.
        """
        back_transcribed_seq = fhe.zeros(len(self))
        for i in range(len(self)):
            back_transcribed_seq[i] = _SeqAbstractBaseClass._backTranscriptionTable[self._data[i]]
        if inplace:
            if self.__class__ == FheSeq:
                raise TypeError("Sequence is immutable")
            self._data = back_transcribed_seq
            return self
        else:
            return self.__class__(back_transcribed_seq)

    def translate(self, table="Standard"):
        """Turn a nucleotide sequence into a protein sequence by creating a new sequence object.

        This method will translate DNA or RNA sequences. It should not
        be used on protein sequences as any result will be biologically
        meaningless.

        Arguments:
         - table - Which codon table to use?  This can be either a name
           (string), or a fhe.LookupTable object (useful for non-standard genetic codes).
           This defaults to the "Standard" table.
        """
        n = len(self)

        if n % 3 != 0:
            raise Error(
                "Partial codon, len(sequence) not a multiple of three. "
            )

        translationTable = fhe.LookupTable(SeqWrapper.get_translationTable(table))

        # first of all, reduce the integers of letters 'ACGU' (and T treated as a U) to 0,1,2,3
        reduced_integers = fhe.zeros(len(self))
        for i in range(len(self)):
            reduced_integers[i] = _SeqAbstractBaseClass._translationReductionTable[self._data[i]]

        protein_seq = fhe.zeros(n//3)

        for i in range(n//3):
            # compute codon index from first, second and third letters
            codon_index = reduced_integers[i*3]*16 + reduced_integers[i*3+1]*4 + reduced_integers[i*3+2]
            protein_seq[i] = translationTable[codon_index]

        return self.__class__(protein_seq)


class FheSeq(_SeqAbstractBaseClass):

    def __init__(self, data, length=None):
        super().__init__(data, length)


class FheMutableSeq(_SeqAbstractBaseClass):

    def __init__(self, data, length=None):
        super().__init__(data, length)

    def reverse(self):
        s=len(self)
        reverse_seq = fhe.zeros(s)
        for i in range(s):
            reverse_seq[i] = self._data[s-i-1]
        self._data = reverse_seq
