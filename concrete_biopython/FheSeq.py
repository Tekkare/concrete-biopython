from abc import ABC
from concrete import fhe
import numpy as np
import numbers

from concrete_biopython.SeqWrapper import SeqWrapper


class _FheSeqAbstractBaseClass(ABC):
    """
    The FHE version Bio.Seq._SeqAbstractBaseClass abstract class

    Methods from Bio.Seq._SeqAbstractBaseClass that cannot be implemented in fhe:
    split, rsplit, strip, lstrip, rstrip


    Dev notes :
    -----------
    
    TODO:

    # Hard to implement in fhe (if doable):
    count, count_overlap , __contains__ , find , rfind , index , rindex, replace

    # Cannot be implemented, because it collides with array __add__ function:
    __radd__

    # Require compatible implementaion of numpy.repeat for concrete to be done properly
    __mul__ , __rmul__ , __imul__

    # Unused for now because all characters are uppercase ( this lowers the letter variable to 5 bits only)
    #   but lower case characters could easily be added to SeqWrapper.LETTERS
    upper , lower , isupper , islower

    # If the class is extended to include undefined sequences:
    defined , defined_ranges


    Possible optimization:
    Child class FheBaseSeq that only encode DNA and RNA bases with integers in range 0..3, which will take only 2 bits for quicker processing

    """

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
                #self._data = _UndefinedSequenceData(length)
                raise NotImplemented
        elif isinstance(data, fhe.tracing.tracer.Tracer):
            if data.size>1:
                self._data = data[:] # take a copy
            elif data.size==1:
                try:
                    self._data = data[:].reshape(1) # take a copy if array of size 1
                except:
                    self._data = data.reshape(1) # just copy value if isolated value
            else:
                raise ValueError("data length must be positive")
        elif isinstance(data, _FheSeqAbstractBaseClass):
            if data._data.size>1:
                self._data = data._data[:] # take a copy
            elif data._data.size==1:
                self._data = data._data[:].reshape(1)  # take a copy
            else:
                self._data = fhe.zeros(0)
        else:
            raise TypeError(
                "data should be a concrete.fhe.tracing.tracer.Tracer object"
            )    


    def toArray(self):
        return self._data

    def __eq__(self, other):
        """Compare the sequence to another sequence or an array.

        Sequences are equal to each other if their sequence contents is identical
        """
        if isinstance(other, _FheSeqAbstractBaseClass):
            other = other._data
        elif isinstance(other, fhe.tracing.tracer.Tracer):
            pass
        else:
            return NotImplemented

        if len(self) != other.size:
            return 0 # return False if the lengths are different
        else:
            return (len(self) - np.sum(self._data == other))==0


    def __lt__(self, other):
        """ Computes wether a sequence is lower than another, in alphabetical order

        A sequence A is lower than a sequence B if and only if:
        there exist an index i such that: A[i] < B[i]  AND  for all k<i, A[k] <= B[k]

        This is the mathematical converse of A >= B (see self.__ge__) which is easier to do in fhe

        If sequence have different length, the shorter one is considered to have extra empty characters
        where the empty character is the lowest in alphabetical order
        """
        return 1-(self >= other)

    def __le__(self, other):
        """ Computes wether a sequence is lower or equal than another, in alphabetical order

        The proposition "A <= B" is equivalent to "B >= A" (see self.__ge__)

        If sequence have different length, the shorter one is considered to have extra empty characters
        where the empty character is the lowest in alphabetical order
        """       
        return other >= self

    def __gt__(self, other):
        """ Computes wether a sequence is greater than another, in alphabetical order

        The proposition "A > B" is the mathematical converse of "B >= A" (see self.__ge__)

        If sequence have different length, the shorter one is considered to have extra empty characters
        where the empty character is the lowest in alphabetical order
        """
        return 1-(other >= self)

    def __ge__(self, other):
        """ Computes wether a sequence is greater or equal than another, in alphabetical order

        A sequence A is greater or equal than a sequence B if and only if:
        for all index i:  either  A[i] >= B[i]  or  there is an index k<i where A[k] > B[k]

        If sequence have different length, the shorter one is considered to have extra empty characters
        where the empty character is the lowest in alphabetical order
        """       

        A=self._data[:]  # copy array
        B=other._data[:] # copy array

        # prepare arrays if they have different size
        if A.size != B.size:
            # append a zero to the shortest sequence and crop the other one
            #   zero is always <= to any other character so no need to append more zeros
            diff= B.size - A.size
            if diff > 0:
                A = np.concatenate((A, fhe.zeros(1).reshape(1)), axis=0)
                B = B[0:A.size]
            else:
                B = np.concatenate((B, fhe.zeros(1).reshape(1)), axis=0)
                A = A[0:B.size]
        n = A.size

        if n==0:
            return 1 # special case if both arrays are empty

        ## compute A[i] >= B[i] array
        ge_array = A >= B
        
        ## compute A[i] > B[i] array
        gt_array = A > B 

        ## compute "there is an index k<i where A[i] > B[i]" array in cum_gt_array:
        # if gt_array[k] is 1 for some k, cum_gt_array[i] will be 1 for all i>k
        cum_gt_array = fhe.zeros(n)
        if n >1:
            cum_gt_array[1] = gt_array[0]
        else:
            return ge_array[0] # special case if array has size one

        for i in range(2,n):
            cum_gt_array[i] = cum_gt_array[i-1] | gt_array[i-1]

        ## now compute " A[i] >= B[i]  or  there is an index k<i where A[i] > B[i] " array in or_array:
        or_array = ge_array | cum_gt_array

        ## return wether or_array is true for all indices
        return (n - np.sum(or_array))==0


    """
    Dev:
    A more elegant alternative to __ge__ function, but way slower :

    def _arr_ge(arr1, arr2, i):
        if (arr2.size - i)== 0:
            return 1
        elif (arr1.size - i) == 0:
            return 0
        else:
            return (arr1[i]>arr2[i]) | ((arr1[i]==arr2[i]) & _FheSeqAbstractBaseClass._arr_ge(arr1, arr2, i+1))

    def __ge__(self, other):
        return _FheSeqAbstractBaseClass._arr_ge(self._data, other._data, 0)

    """

    """
    Dev:
    A second alternative to __ge__ function, quite good but still slower than original:
    
    def __ge__(self, other):

        A=self._data[:]  # copy array
        B=other._data[:] # copy array

        # prepare arrays if they have different size
        if A.size != B.size:
            # append a zero to the shortest sequence and crop the other one
            #   zero is always <= to any other character so no need to append more zeros
            diff= B.size - A.size
            if diff > 0:
                A = np.concatenate((A, fhe.zeros(1).reshape(1)), axis=0)
                B = B[0:A.size]
            else:
                B = np.concatenate((B, fhe.zeros(1).reshape(1)), axis=0)
                A = A[0:B.size]
        n = A.size

        arr_g = A>B
        arr_e = A==B

        ge = 1
        for i in range(n-1,-1,-1):
            ge = arr_g[i] | (arr_e[i] & ge)

        return ge

    """    

    def __len__(self):
        """Return the length of the sequence."""
        return self._data.size

    def __iter__(self):
        """Return an iterable of the sequence."""
        return iter(self._data)

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
        if isinstance(other, _FheSeqAbstractBaseClass):
            return self.__class__(np.concatenate((self._data, other._data), axis=0)) # np.concatenate will deep copy arrays
        elif isinstance(other, fhe.tracing.tracer.Tracer):
            return self.__class__(np.concatenate((self._data, other), axis=0))
        else:
            return NotImplemented

    def startswith(self, prefix, start=None, end=None):
        """Return True if data starts with the specified prefix, False otherwise.

        With optional start, test data beginning at that position.
        With optional end, stop comparing data at that position.
        prefix can also be a tuple of bytes to try.
        """
        if isinstance(prefix,_FheSeqAbstractBaseClass):
            l=len(prefix)
        elif isinstance(prefix, fhe.tracing.tracer.Tracer):
            l=prefix.size
        else:
            return NotImplemented

        if start and end:
            return self[start:end] == prefix
        elif start:
            return self[start:start+l] == prefix
        elif end:
            return self[0:np.min((l,end))] == prefix
        else:
            return self[0:l] == prefix

    def endswith(self, suffix, start=None, end=None):
        """Return True if data ends with the specified suffix, False otherwise.

        With optional start, test data beginning at that position.
        With optional end, stop comparing data at that position.
        suffix can also be a tuple of bytes to try.
        """
        if isinstance(suffix,_FheSeqAbstractBaseClass):
            l=len(suffix)
        elif isinstance(suffix, fhe.tracing.tracer.Tracer):
            l=suffix.size
        else:
            return NotImplemented

        if start and end:
            return self[start:end] == suffix
        elif start:
            return self[np.max((start,len(self)-l)):] == suffix
        elif end:
            return self[end-l:end] == suffix
        else:
            return self[len(self)-l:] == suffix

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
        reduced_integers = _FheSeqAbstractBaseClass._translationReductionTable[self._data]

        protein_seq = fhe.zeros(n//3)

        for i in range(n//3):
            # compute codon index from first, second and third letters
            codon_index = reduced_integers[i*3]*16 + reduced_integers[i*3+1]*4 + reduced_integers[i*3+2]
            protein_seq[i] = translationTable[codon_index]

        return self.__class__(protein_seq)


    def complement(self, inplace=False):
        """Return the complement as an RNA sequence.

        Any T in the sequence is treated as a U
        Any other letter than ACGT will be changed to the character with index 0 which should not be done

        The sequence is modified in-place and returned if inplace is True

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if
        ``complement_rna`` is called on a ``Seq`` object with ``inplace=True``.
        """
        complement_seq = _FheSeqAbstractBaseClass._DNAcomplementTable[self._data]
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
        complement_seq = _FheSeqAbstractBaseClass._RNAcomplementTable[self._data]
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
        transcribed_seq = _FheSeqAbstractBaseClass._transcriptionTable[self._data]
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
        back_transcribed_seq = _FheSeqAbstractBaseClass._backTranscriptionTable[self._data]
        if inplace:
            if self.__class__ == FheSeq:
                raise TypeError("Sequence is immutable")
            self._data = back_transcribed_seq
            return self
        else:
            return self.__class__(back_transcribed_seq)

    def join(self, other):
        """Return a merge of the sequences in other, spaced by the sequence from self.

        Accepts a FheSeq object, FheMutableSeq object, or array (and iterates over
        the integers), or an iterable containing Seq, MutableSeq, or string
        objects. These arguments will be concatenated with the calling sequence
        as the spacer:
        """
        joindata = []
        if isinstance(other, _FheSeqAbstractBaseClass):
            for i in range(len(other)):
                joindata.append(other[i].reshape(1))
        elif isinstance(other, fhe.tracing.tracer.Tracer):
            for i in range(other.size):
                joindata.append(other[i].reshape(1))
        elif isinstance(other, list) or isinstance(other, tuple):
            for data in other:
                if isinstance(data, _FheSeqAbstractBaseClass):
                    joindata.append(data.toArray())
                elif isinstance(data, fhe.tracing.tracer.Tracer):
                    joindata.append(data)
                else:
                    raise NotImplemented
        else:
            NotImplemented

        concatenation = joindata[0]
        for i in range(1,len(joindata)):
            concatenation = np.concatenate((concatenation,self._data,joindata[i]),axis=0)

        return self.__class__(concatenation)


class FheSeq(_FheSeqAbstractBaseClass):
    """
    The FHE version Bio.Seq.Seq class

    Dev notes:
    ----------
    
    TODO:

    __hash__

    """
    def __init__(self, data, length=None):
        super().__init__(data, length)


class FheMutableSeq(_FheSeqAbstractBaseClass):
    """
    The FHE version Bio.Seq.MutableSeq class

    This method from Bio.Seq.FheMutableSeq that cannot be implemented in fhe: remove
    """
    def __init__(self, data, length=None):
        super().__init__(data, length)

    def reverse(self):
        s=len(self)
        reverse_seq = fhe.zeros(s)
        for i in range(s):
            reverse_seq[i] = self._data[s-i-1]
        self._data = reverse_seq

    def __setitem__(self, index, value):
        """Set a subsequence of single letter via value parameter.
        """
        if isinstance(index, numbers.Integral):
            # Replacing a single letter with a new string
            self._data[index] = value
        else:
            # Replacing a sub-sequence
            if isinstance(value, _FheSeqAbstractBaseClass):
                self._data[index] = value._data
            elif isinstance(value, fhe.tracing.tracer.Tracer):
                self._data[index] = value
            else:
                raise TypeError(f"received unexpected type '{type(value).__name__}'")

    def __delitem__(self, index):
        """Delete a subsequence of single letter.
        """
        if isinstance(index, numbers.Integral):
            # Replacing a single letter with a new string
            if(index > 0 or index <-1):
                self._data = (self[0:index]+self[index+1:])._data
            elif index==0:
                self._data = self[1:]._data
            else: #-1
                self._data = self[0:index]._data
        else:
            self._data = (self[0:index.start]+self[index.stop:])._data

    def append(self, c):
        """Add a single letter to the mutable sequence object.
        No return value.
        """
        if not isinstance(c, fhe.tracing.tracer.Tracer):
            raise TypeError('c must be of type fhe.tracing.tracer.Tracer')
        if c.size > 1:
            raise ValueError('c must be of size 1')
        self._data = (self+c.reshape(1))._data

    def insert(self, i, c):
        """Add a single letter to the mutable sequence object at a given index.
        No return value.
        """
        if not isinstance(c, fhe.tracing.tracer.Tracer):
            raise TypeError('c must be of type fhe.tracing.tracer.Tracer')
        if c.size > 1:
            raise ValueError('c must be of size 1')
        self._data = (self[0:i]+c.reshape(1)+self[i:])._data

    def pop(self, i=-1):
        """Remove a subsequence of a single letter at given index.
        Returns the last character of the sequence.
        """
        c = self._data[i]
        self.__delitem__(i)
        return c
    
    def extend(self, other):
        """Add a sequence to the original mutable sequence object.
        No return value.
        """
        if isinstance(other, _FheSeqAbstractBaseClass) or isinstance(other,fhe.tracing.tracer.Tracer):
            self._data = (self+other)._data
        else:
            raise TypeError("expected a string, Seq or MutableSeq")

