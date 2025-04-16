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
            
def kmers (sequence, k):
    kmers = {}
    for pos in range(len(sequence)-k+1):
        kmer = sequence[pos:pos+k]
        if not kmer in kmers:
            kmers[kmer] = 1
        else:
            kmers[kmer] += 1
    return kmers    

def graphe_de_de_Bruijn(filename):
    for read in read_gz(filename):
        read_seq = str(read.seq)
        list_kmers = kmers(read_seq, 3)
    pass


if __name__ == '__main__':
    s = "AAAA"
    t = "ACCC"
    suff = s[-1]
    if t.startswith(suff):
        print ("oui")
    else : 
        print ("non")

    print(kmers("AAAAATTCCCCTTCTGGG", 5))