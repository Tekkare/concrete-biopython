import numpy as np
from Bio.Seq import Seq, MutableSeq


class SeqWrapper():

    # create type letter related variables
    LETTERS = '*ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    LETTERS_TO_INTEGERS = {letter: index for index, letter in enumerate(LETTERS)}
    INTEGERS_TO_LETTERS = {index: letter for index, letter in enumerate(LETTERS)}

    def __init__(self, data, length=None):
        """
        Arguments:
        - data: can be a BioPython Seq, MutableSeq, str, list or numpy array
        """

        # if data is an integer array, translate it to letters if needed
        if isinstance(data, np.ndarray) or isinstance(data, list):
            valid_integers = set(SeqWrapper.INTEGERS_TO_LETTERS.keys())
            if not set(data).issubset(valid_integers):
                raise ValueError("Invalid integers, must be in range 0:",len(SeqWrapper.LETTERS))
            data = ''.join([SeqWrapper.INTEGERS_TO_LETTERS[integer] for integer in data])

        elif isinstance(data, Seq) or isinstance(data, MutableSeq):
            data = data.__str__()

        # Verify sequence letters
        valid_chars = set(SeqWrapper.LETTERS_TO_INTEGERS.keys())
        if not set(data).issubset(valid_chars):
            raise ValueError("Invalid characters, must be in ", SeqWrapper.LETTERS)

        self._data = data

    def __str__(self):
        return self._data

    def toSeq(self):
        return Seq(self._data)

    def toMutableSeq(self):
        return MutableSeq(self._data)        

    def toIntegers(self):
        # Map letters to integers
        return np.array([SeqWrapper.LETTERS_TO_INTEGERS[char] for char in self._data])

    def getLetters():
        return str(SeqWrapper.LETTERS)

    def _makeTable(letter_mapping):
        """
        Makes a integer mapping table from a letter mapping
        Absent letters are mapped to the index 0 which should not be done
        """
        integer_mapping = { SeqWrapper.LETTERS_TO_INTEGERS[letter]: SeqWrapper.LETTERS_TO_INTEGERS[letter_mapping[letter]] for letter in letter_mapping.keys() }
        return [ integer_mapping[i] if i in integer_mapping else 0 for i in range(len(SeqWrapper.LETTERS)) ]

    def get_DNA_complementTable():
        """
        Computes an integer table to map DNA letters as integers to their complement letter as integer
        Any U is treated like a T
        Other letters are mapped to the integer 0, which should not be done
        """
        return SeqWrapper._makeTable( {'A':'T', 'T':'A', 'U':'A', 'C':'G', 'G':'C'} )

    def get_RNA_complementTable():
        """
        Computes an integer table to map RNA letters as integers to their complement letter as integer
        Any T is treated like a U
        Other letters are mapped to the integer 0, which should not be done
        """
        return SeqWrapper._makeTable( {'A':'U', 'T':'A', 'U':'A', 'C':'G', 'G':'C'} )

    def get_transcriptionTable():
        """
        Computes an integer table to map DNA letters as integers to their transcribed RNA letter as integer
        Only T is changed to U, ACG are unchanged
        Any U is treated like a T
        Other letters are mapped to the integer 0, which should not be done
        """
        return SeqWrapper._makeTable( {'A':'A', 'T':'U', 'U':'U', 'C':'C', 'G':'G'} )

    def get_back_transcriptionTable():
        """
        Computes an integer table to map RNA letters as integers to their transcribed DNA letter as integer
        Only U is changed to T, ACG are unchanged
        Any T is treated like a U
        Other letters are mapped to the integer 0, which should not be done
        """
        return SeqWrapper._makeTable( {'A':'A', 'U':'T', 'T':'T', 'C':'C', 'G':'G'} )

    def get_translationReductionTable():
        """
        Computes a table to map the integer version of RNA letters 'ACGU' to intgers in 0..3
            which is required to process codons (groups of letters)
        Letter T is treated like a U
        Other letters would be mapped to 0, which should not be done
        """
        return SeqWrapper._makeTable( {'A':SeqWrapper.LETTERS[0], 'U':SeqWrapper.LETTERS[3], 'T':SeqWrapper.LETTERS[3], 'C':SeqWrapper.LETTERS[1], 'G':SeqWrapper.LETTERS[2]} )

    def get_translationTable(table="Standard"):
        """
        Computes a table to map RNA letters as integers in 0..3 to their translated protein letter as integer
        This default table corresponds to the standard translation table of DNA
        Any T is treated like a U
        Other letters are mapped to the integer 0, which should not be done
        """
        import itertools
        codons = [''.join(comb) for comb in itertools.product('ACGU', repeat=3)]
        codons_translated = [Seq(c).translate(table).__str__() for c in codons]

        return [SeqWrapper.LETTERS_TO_INTEGERS[letter] for letter in codons_translated]

