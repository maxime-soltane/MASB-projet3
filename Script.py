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

def write_fasta (seq, filename = "assemble_seq.fasta"):
    with open (filename, "w") as f :
        f.write("> assembled sequence\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+60] + "\n")


if __name__ == '__main__':
    # Test with the longer sequence
    f = read_gz("Level0.fa.gz")
    kmers_dict = {}
    for seq in f:
        kmers_dict.update(kmers(str(seq.seq), 21))

    longer_sequence = "CCCACGGACGCCAGAACGGGCGTTCTCCCTAGCGTGCGCCCTGCAGAACGTTCGCGAGAACGACAGAACTCACGGACGTTCTCCCTATCGACCGTGCGCAAGAACGTCCGGCCGTACGCCCTATAGAACGAGCGCCCGCTCGGCCGTGCTATAGAACTCTCGGCCTCACGGAAGAACGTTCGTGCGCCAGAACTATCTCACGCCCTAAAGTG"
    
    dbg = DeBruijnGraph(kmers_dict)
    assembled_long_seq = dbg.assemble_sequence()
    print("\nLong sequence test:")
    print(f"Original length: {len(longer_sequence)}")
    print(f"Assembled length: {len(assembled_long_seq)}")
    print(f"Correct assembly: {longer_sequence == assembled_long_seq}")

    #Test with level1
    print("")
    print("Test with level1")
    start = time()
    f2 = read_gz("Level1.fa.gz")

    kmers_dict2 = defaultdict(int)  
    for seq in f2:
        for kmer, count in kmers(str(seq.seq), 21).items():
            kmers_dict2[kmer] += count  
    print(f"Nombre de kmers avant filtrage : {len(kmers_dict2)}")

    f_kmers = kmers_filter(kmers_dict2, 3)
    print(f"Nombre de kmers restant après filtrage : {len(f_kmers)}")

    dbg2 = DeBruijnGraph(f_kmers)
    assembled_seq2 = dbg2.assemble_sequence()
    print(f"250 premiers nucléotides de la séquence assemblé : {assembled_seq2[0:250]}")
    print(f"Taille de la séquence assemblée : {len(assembled_seq2)}")
    
    end = time()
    print(f"Temps d'exécution sur level1 : {end-start}")

    write_fasta(assembled_seq2)