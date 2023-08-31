import numpy as np
import itertools
from Bio.Seq import Seq, MutableSeq


class SeqWrapper():
    """
    A wrapper class to interface Concrete-BioPython with BioPython

    Dev notes:
    ---------
    New letters can be added to the default alphabet, and tables will adapt automatically.
    LETTERS will be automatically sorted alphabetically, so that the > operator corresponds to the one in Bio.Seq.
    Adding more characters will increase the number of bits required to encode for each character.
    The empty character character \0 should always be kept in the default alphabet
    """

    # All letters of the default alphabet
    _DEFAULT_ALPHABET = '\0*ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    # Initialize letter related variables 
    LETTERS = ''.join((sorted(set(_DEFAULT_ALPHABET))))
    LETTERS_TO_INTEGERS = {letter: index for index, letter in enumerate(LETTERS)}
    INTEGERS_TO_LETTERS = {index: letter for index, letter in enumerate(LETTERS)}

    def __init__(self, *args): 
        raise Exception('This class should not be instanciated')

    def setAlphabet(alphabet):
        """
        Set a custom alphabet for optimized computations, typically for DNA/RNA only
        """
        if not isinstance(alphabet, str):
            raise ValueError('alphabet must be a string')
        if not '\0' in alphabet:
            alphabet += '\0'
        SeqWrapper.LETTERS = ''.join((sorted(set(alphabet))))
        SeqWrapper.LETTERS_TO_INTEGERS = {letter: index for index, letter in enumerate(SeqWrapper.LETTERS)}
        SeqWrapper.INTEGERS_TO_LETTERS = {index: letter for index, letter in enumerate(SeqWrapper.LETTERS)}

    def resetAlphabet():
        """
        Reset alphabet to default value
        """
        SeqWrapper.setAlphabet(SeqWrapper._DEFAULT_ALPHABET)

    def toStr(array):
        """
        Map integer array to string
        """
        # array is an integer array, translate it to letters 
        if isinstance(array, np.ndarray):
            valid_integers = set(SeqWrapper.INTEGERS_TO_LETTERS.keys())
            if not set(array).issubset(valid_integers):
                raise ValueError("Invalid integers, must be in range 0:",len(SeqWrapper.LETTERS))
            string = ''.join([SeqWrapper.INTEGERS_TO_LETTERS[integer] for integer in array])

        else:
            raise ValueError('input must be an integer array')

        return string

    def toSeq(array):
        """
        Map integer array to Seq object
        """        
        return Seq(SeqWrapper.toStr(array))

    def toMutableSeq(array):
        """
        Map integer array to MutableSeq object
        """        
        return MutableSeq(SeqWrapper.toStr(array))

    def toIntegers(data):
        """
        Maps letters to integers

        Arguments:
        - data: can be a BioPython Seq, MutableSeq, or str
        """
        if isinstance(data, Seq) or isinstance(data, MutableSeq):
            data = data.__str__()

        elif isinstance(data, str):
            data = data

        else:
            raise ValueError('data must be either BioPython Seq, MutableSeq or str ')

        # Verify sequence letters
        valid_chars = set(SeqWrapper.LETTERS_TO_INTEGERS.keys())
        if not set(data).issubset(valid_chars):
            raise ValueError("Invalid characters, must be in ", SeqWrapper.LETTERS)

        # return array
        return np.array([SeqWrapper.LETTERS_TO_INTEGERS[char] for char in data])

    def maxInteger():
        """
        Returns the maximum integer, encoding the last character from the alphabet LETTERS
        """        
        return len(SeqWrapper.LETTERS)-1

    def _makeTable(letter_mapping):
        """
        Makes an integer mapping table from a letter mapping
        Absent letters are mapped to the index 0 which should not be done
        """
        try:
            integer_mapping = { SeqWrapper.LETTERS_TO_INTEGERS[letter]: SeqWrapper.LETTERS_TO_INTEGERS[letter_mapping[letter]] for letter in letter_mapping.keys() }
        except KeyError as e: # if the alphabet was changed and does not contain the right letters anymore
            raise KeyError('SeqWrapper alphabet has been set without the required letter:' + str(e))

        return [ integer_mapping[i] if i in integer_mapping else 0 for i in range(len(SeqWrapper.LETTERS)) ]

    def get_DNA_complementTable():
        """
        Computes an integer table to map DNA letters (as integers) to their complement letter (as integer)
        Any U is treated like a T
        Other letters are mapped to the integer 0, which should not be done
        """
        return SeqWrapper._makeTable( {'A':'T', 'T':'A', 'U':'A', 'C':'G', 'G':'C'} )

    def get_RNA_complementTable():
        """
        Computes an integer table to map RNA letters (as integers) to their complement letter (as integer)
        Any T is treated like a U
        Other letters are mapped to the integer 0, which should not be done
        """
        return SeqWrapper._makeTable( {'A':'U', 'T':'A', 'U':'A', 'C':'G', 'G':'C'} )

    def get_transcriptionTable():
        """
        Computes an integer table to map DNA letters (as integers) to their transcribed RNA letter (as integer)
        Only T is changed to U, ACG are unchanged
        Any U is treated like a T
        Other letters are mapped to the integer 0, which should not be done
        """
        return SeqWrapper._makeTable( {'A':'A', 'T':'U', 'U':'U', 'C':'C', 'G':'G'} )

    def get_back_transcriptionTable():
        """
        Computes an integer table to map RNA letters (as integers) to their transcribed DNA letter (as integer)
        Only U is changed to T, ACG are unchanged
        Any T is treated like a U
        Other letters are mapped to the integer 0, which should not be done
        """
        return SeqWrapper._makeTable( {'A':'A', 'U':'T', 'T':'T', 'C':'C', 'G':'G'} )

    def get_translationReductionTable():
        """
        Computes a table to map the integer version of RNA letters 'ACGU' to integers in 0..3
            which is required to process codons (groups of letters)
        Letter T is treated like a U
        Other letters would be mapped to 0, which should not be done
        """
        return SeqWrapper._makeTable( {'A':SeqWrapper.LETTERS[0], 'U':SeqWrapper.LETTERS[3], 'T':SeqWrapper.LETTERS[3], 'C':SeqWrapper.LETTERS[1], 'G':SeqWrapper.LETTERS[2]} )

    def get_translationTable(table="Standard"):
        """
        Computes a table to map RNA letters (as integers) in 0..3 to their translated protein letter (as integer)
        This default table corresponds to the standard translation table of DNA
        Any T is treated like a U
        Other letters are mapped to the integer 0, which should not be done
        """
        # generate all possible codons in alphabetical order
        codons = [''.join(comb) for comb in itertools.product('ACGU', repeat=3)]

        # translate the codons into a protein letter with BioPython Seq object
        codons_translated = [Seq(c).translate(table).__str__() for c in codons]

        # return a table of the integer represntation of the protein letters
        try:
            return [SeqWrapper.LETTERS_TO_INTEGERS[letter] for letter in codons_translated]    
        except KeyError as e: # if the alphabet was changed and does not contain the right letters anymore
            raise KeyError('SeqWrapper alphabet has been set without the required letter:' + str(e))
        

