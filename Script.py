import gzip
from Bio import SeqIO   
from DeBruijnGraph import * 
    
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
            
def kmers (sequence, k):
    kmers = {}
    for pos in range(len(sequence)-k+1):
        kmer = sequence[pos:pos+k]
        if not kmer in kmers:
            kmers[kmer] = 1
        else:
            kmers[kmer] += 1
    return kmers   

def kmers_filter (kmers_dict, threshold = 1):
    """
    
    >>> k_dict = {"A" : 4, "B" : 1, "C" : 12, "D" : 5}
    >>> kmers_filter(k_dict, 5)
    {'C': 12, 'D': 5}
    """
    f_kmers = {}
    for kmer, count in kmers_dict.items():
        if count >= threshold:
            f_kmers[kmer] = count
    
    return f_kmers

if __name__ == '__main__':
    # Test with a simple sequence
    test_sequence = "ATGGGTGGTGGTATG"
    kmers_dict = kmers(test_sequence, 3)

    dbg = DeBruijnGraph(kmers_dict)
    print(f"Graph: {dict(dbg.get_graph())}")
    print(f"Start node: {dbg.get_start()}")
    print(f"Simple path: {dbg.simple_path()}")
    assembled_seq = dbg.assemble_sequence()
    print(f"Assembled sequence: {assembled_seq}")
    print(f"Original sequence: {test_sequence}")
    print(f"Correct assembly: {test_sequence == assembled_seq}")

    #Ne fonctionne pas a-t-on des répétitions dans les séquences fournies

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
    f2 = read_gz("Level1.fa.gz")

    kmers_dict2 = {}
    for seq in f2:
        kmers_dict2.update(kmers(str(seq.seq), 21))
    print(len(kmers_dict2))
    f_kmers = kmers_filter(kmers_dict2)
    print(len(f_kmers))

    dbg2 = DeBruijnGraph(f_kmers)
    assembled_seq2 = dbg2.assemble_sequence()
    print(assembled_seq2)
    #Une séquence est bien assemblée voir si c'est la bonne mais tous les kmers ne sont vu qu'une fois 