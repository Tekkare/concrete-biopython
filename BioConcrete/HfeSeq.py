from abc import ABC
from concrete import fhe
import numpy as np

class _SeqAbstractBaseClass(ABC):

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

    def complement(self, inplace=False):
        complement_seq = fhe.zeros(self._data.size)
        for i in range(self._data.size):
            complement_seq[i] = 3-self._data[i]
        if inplace:
            self._data = complement_seq
            return self
        else:
            return self.__class__(complement_seq)

    def reverse(self):
        s=self._data.size
        reverse_seq = fhe.zeros(s)
        for i in range(s):
            reverse_seq[i] = self._data[s-i-1]
        self._data = reverse_seq

    def reverse_complement(self, inplace=False):
        reverse_complement = self.complement(inplace)
        reverse_complement.reverse()
        return reverse_complement


class HfeSeq(_SeqAbstractBaseClass):

    def __init__(self, data, length=None):
        super().__init__(data, length)


class HfeMutableSeq(_SeqAbstractBaseClass):

    def __init__(self, data, length=None):
        super().__init__(data, length)
