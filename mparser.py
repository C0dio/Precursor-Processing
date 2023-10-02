import gzip
from bio import Bio

class File:
    """
    File Class

    simple class to showcase alternative methods to reading file contents &
    the process of compressing and decompressing gzipped files.
    """
    @staticmethod
    def convert_to_gzip(target: str, output: str):
        contents = None
        with open(target, 'rb') as tfile:
            contents = tfile.read()

        with gzip.open(output, 'wb') as ofile:
            ofile.write(contents)

    @staticmethod
    def convert_to_tsv(target: str, output: str):
        # get the contents from the file
        contents = None
        with gzip.open(target, "rt") as tfile:
            contents = tfile.read().split('\n')

        # convert to iterator & grab the first item
        itercontents = iter(contents)
        id = next(itercontents).replace('>', '')
        
        # traverse the file and write to output
        with open(output, 'wt') as ofile:
            while id:
                # get the seq
                sequence = next(itercontents)

                # get the notation and MFE
                notation_mfe = next(itercontents).split(' ')
                # format the MFE without brackets
                notation_mfe[1] = notation_mfe[1].strip('()')

                # write to the file
                ofile.write(f'{id}\t{sequence}\t{notation_mfe[0]}\t{notation_mfe[1]}\n')

                # get the next id
                id = next(itercontents).replace('>', '')


class Feature:
    """
    Feature Class

    used to store information about a feature, typically obtained through gff3 file, and
    contains consistent functionality and a method stub for classes to override.  
    """
    headers = {}

    def __init__(self, seqid: str, start: str, end: str, transcript: str):
        self.seqid: str = seqid
        self.start: int = int(start)
        self.end: int = int(end)
        self.transcript: str = transcript

    def precursor_sequence(self, chromosomes) -> str:
        # get the correct chromosome, default is index 7
        chromosome = chromosomes.get(self.seqid, '')
        # skip the first line about chromosome details
        chromosome = "".join(chromosome.split("\n")[1:])
        # return the sequence between the start-end indexes
        return self.extract_chr(chromosome)

    def extract_chr(self, _):
        raise NotImplementedError()

class FFeature(Feature):
    "Forward Feature Class"
    def extract_chr(self, chromosome):
        return chromosome[self.start - 1: self.end]

class BFeature(Feature):
    "Backward Feature Class"
    def extract_chr(self, chromosome):
        return Bio.revcomp(chromosome[self.start - 1: self.end])


def convert_text_to_feature(line: str) -> Feature:
    """
    Read each line of the gff3 file and convert them into a Feature instance
    """
    columns = line.split("\t")
    if columns[6] == '+':
        return FFeature(
            columns[0], # seqid
            columns[3], # start
            columns[4], # end
            columns[8]  # transcript
        )
    return BFeature(
        columns[0],
        columns[3],
        columns[4],
        columns[8]
    )


def convert_file_to_features(path: str) -> [Feature]:
    """
    Read a gff3 file and convert its contents into a Feature list
    """
    features: [Feature] = []
    # read the gz file as a string and split into lines
    with gzip.open(path, 'rt') as file:
        lines = file.read().split("\n")
        # read through the file one line at a time
        for line in lines:
            if "pre_miRNA" in line:
                # convert the line to a feature & add to list
                features.append(convert_text_to_feature(line))
    return features


def nsv_to_tsv():
    """
    simplified example demonstrating the compression to and from a .gz file
    """
    File.convert_to_gzip('rna.nsv', 'rna.gz')
    File.convert_to_tsv('rna.gz', 'rna.tsv')


def main(compress = False):
    if compress:
        File.convert_to_gzip('inputs/input.gff3', 'inputs/input.gff3.gz')
        File.convert_to_gzip('inputs/input.fa', 'inputs/input.fa.gz')

    # read in the file
    # NB: such files can be obtained though ensemble's api or website
    features = convert_file_to_features('inputs/input.gff3.gz')

    # read the FASTA file
    results = []
    with gzip.open('inputs/input.fa.gz', 'rt') as fasta:
        fasta = fasta.read()

        # split the FASTA file into its respective strands/chromosomes
        chromosomes = {}
        for chromosome in fasta.split(">")[1:]:
            chromosomes.update(Bio.format_fasta_chromosome(chromosome)) 

        # read the precursor sequence based on the type of Feature
        for feature in features:
            # compile all the information and save it
            # TODO: add padding function to make id consistently padded from precursor seq.
            results.append(
                f"{feature.seqid}\t\
                  {feature.precursor_sequence(chromosomes)}"
            )

    # do something with the results
    for result in results:
        print(result)
    
    # print some simple verification
    print("".ljust(150, '-'))
    print(f"Expected: {len(features)} sequences, Actual: {len(results)}. Success: {len(features) == len(results)}")


if __name__ == '__main__':
    main()