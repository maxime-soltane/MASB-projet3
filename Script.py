from pybloom_live import ScalableBloomFilter   
import gzip
from Bio import SeqIO    
    
def read_gz(filename):
        """
        Open a Fasta or FastQ file, compressed (.gz) or not, in reading mode and return sequence's iterator.

        Parameters:
        :param filename: the name of the file to work on

        Raises:
        ValueError: if file format isn't supported
        """
        #Define the correct open function wether the file is compressed or not
        open_func = gzip.open if filename.endswith(".gz") else open  

        with open_func(filename, "rt") as file:
            # Open FastQ file
            if filename.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):  
                yield from SeqIO.parse(file, 'fastq')
            # Open FastA file    
            elif filename.endswith(('.fasta', '.fna', '.fa', '.fasta.gz', '.fna.gz', '.fa.gz')):
                yield from SeqIO.parse(file, 'fasta')
            else:  
                raise ValueError(f"file format not supported: '{filename}'")
            
def bloom_filter(sequence, k) -> ScalableBloomFilter:
    """ 
    Initializes a Bloom filter with all k-mers in the genome

    Returns:
    The Bloom filter containing the genome's kmers
    """
    #Initialize the bloom filter with an estimated capacity and the error rate of the filter returning false positives
    bloom = ScalableBloomFilter(initial_capacity=len(sequence)//k * 2, error_rate=0.01)
        
    # Loop through the sequence to extract and add kmers in the bloom filter
    for i in range(len(sequence) - k + 1):
        bloom.add(sequence[i:i+k])

    return bloom


if __name__ == '__main__':
    k = 'CCC'
    for s in read_gz("Level0.fa.gz"):
        seq = str(s.seq)
        bf = bloom_filter(seq, 3)
        if k in bf:
            print("oui")
        else:
            print("c mor")
