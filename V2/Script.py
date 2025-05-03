#Import des modules nécessaires au fonctionnement du script
import gzip
from Bio import SeqIO   
from collections import defaultdict
from typing import Dict
from matplotlib import pyplot as plt

def read_gz(filename: str):
        """
        Open a Fasta or FastQ file, compressed (.gz) or not, in reading mode and returns sequence's iterator.

        Parameter:
        filename: the name of the file to work on

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
            
def kmers (sequence: str, k: int) -> Dict[str, int]:
    """
    Generates all kmers from a sequence and count their occurrences.

    Parameters:
    sequence: A nucleotidic sequence
    k: The size of the kmer

    Returns:
    A dictionary of kmers and their counts

    Examples:
    >>> kmers("ATCGGCAT", 3)
    {'ATC': 1, 'TCG': 1, 'CGG': 1, 'GGC': 1, 'GCA': 1, 'CAT': 1}
    >>> kmers("AAAAAA", 2)
    {'AA': 5}
    >>> kmers("ATG", 5)
    {}
    """
    kmers = defaultdict(int)
    # Browse every positions where kmers are possibles
    for pos in range(len(sequence)-k+1):
        # Add them to the dictionary with their number of occurrences
        kmers[sequence[pos:pos+k]] += 1
    return kmers

def kmers_filter (kmers_dict: Dict[str, int], threshold: int= 1) -> Dict[str, int]:
    """
    Filters kmers according to a minimum count threshold.

    Parameters:
    kmers_dict: The dictionary of kmer counts
    threshold: The minimum count to keep a kmer

    Returns:
    A filtered dictionary of kmers

    Examples:
    >>> k_dict = {"A" : 4, "B" : 1, "C" : 12, "D" : 5}
    >>> kmers_filter(k_dict, 5)
    {'C': 12, 'D': 5}
    >>> kmers_filter({'A': 2, 'T': 3}, 3)
    {'T': 3}
    """
    f_kmers = defaultdict(int)
    for kmer, count in kmers_dict.items():
        if count >= threshold:
            f_kmers[kmer] = count
    
    return f_kmers

def abundance_hist(kmers_dict: Dict[str, int]) -> None:
    """
    Plots a histogram of kmer abundance on a logarithmic y axis.

    Parameter:
    kmer_dict: The dictionary of kmers and their counts
    """
    abundance_hist = defaultdict(int)
    for count in kmers_dict.values():
            abundance_hist[count] += 1

    x = sorted(abundance_hist.keys())
    y = [abundance_hist[i] for i in x]

    plt.figure(figsize=(12, 6))
    plt.bar(x, y, color="dodgerblue", width=1.0, edgecolor='black')

    plt.title("K-mer abundance histogram")
    plt.xlabel("K-mer multiplicity")
    plt.yscale('log')
    plt.ylabel("Number of Distinct K-mers (log scale)")

    plt.tight_layout()
    plt.show()