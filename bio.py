import gzip

class Bio:
    """
    Bio class used to house simple bio functions because I can't be bothered downloading the package.
    """
    def __init__(self):
        pass

    @staticmethod
    def revcomp(sequence: str) -> str:
        # unwrap the generator object from yield
        return "".join(list(Bio.__revcomp(sequence)))
    
    @staticmethod
    def __revcomp(sequence: str) -> str:
        # return a generator object of the flipped acids
        reverse_seq = sequence[::-1]
        for acid in reverse_seq:
            yield Bio.__flip(acid)
    
    @staticmethod
    def __flip(acid: str) -> str:
        # python 3.10 and higher would have switch/match statements
        # get the compliment of each acid
        if acid == 'A':
            return 'T'
        if acid == 'T':
            return 'A'
        if acid == 'G':
            return 'C'
        return 'G'

    @staticmethod
    def format_fasta_chromosome(text: str) -> (str, str):
        # seperate the header and the chromosome, return as a key value pair
        header, chromosome = tuple(text.split('\n', 1))
        return {header.split(' ')[0]: chromosome}
