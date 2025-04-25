import gzip
from Bio import SeqIO   
from DeBruijnGraph import * 
from collections import defaultdict
from time import time
from typing import Dict
    
def read_gz(filename: str):
        """
        Open a Fasta or FastQ file, compressed (.gz) or not, in reading mode and return sequence's iterator.

        Parameters:
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
    kmers = {}
    for pos in range(len(sequence)-k+1):
        kmer = sequence[pos:pos+k]
        if not kmer in kmers:
            kmers[kmer] = 1
        else:
            kmers[kmer] += 1
    return kmers   

def kmers_filter (kmers_dict: Dict[str, int], threshold: int= 1) -> Dict[str, int]:
    """
    Filters kmers according to a minimum count threshold.

    Parameters:
    kmers_dict: The dictionary of kmer counts
    threshold: The minimum count to keep a kmer

    Returns:
    A filtered dictionnary of kmers

    Examples:
    >>> k_dict = {"A" : 4, "B" : 1, "C" : 12, "D" : 5}
    >>> kmers_filter(k_dict, 5)
    {'C': 12, 'D': 5}
    >>> kmers_filter({'A': 2, 'T': 3}, 3)
    {'T': 3}
    """
    f_kmers = {}
    for kmer, count in kmers_dict.items():
        if count >= threshold:
            f_kmers[kmer] = count
    
    return f_kmers

if __name__ == '__main__':
    #Test with level0
    start = time()
    f = read_gz("Level0.fa.gz")

    kmers_dict = defaultdict(int)
    for seq in f:
        for kmer, count in kmers(str(seq.seq), 21).items():
            kmers_dict[kmer] += count

    dbg = Graph(kmers_dict)
    dbg.get_all_contigs("level0_contig.fa")
    end = time()
    print(f"Temps d'execution sur Level0 = {end-start}\n")

    #Test with level1
    start1 =time()
    f2 = read_gz("Level1.fa.gz")

    kmers_dict2 = defaultdict(int)
    for seq in f2:
        for kmer, count in kmers(str(seq.seq), 31).items():
            kmers_dict2[kmer] += count

    f2_kmers = kmers_filter(kmers_dict2, 5)

    dbg2 = Graph(f2_kmers)
    dbg2.get_all_contigs("level1_contig.fa_21_5")
    end1 = time()
    print(f"Temps d'execution sur Level1 = {end1-start1}\n")