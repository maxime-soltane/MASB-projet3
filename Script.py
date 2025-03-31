import gzip
from Bio import SeqIO    
from collections import defaultdict
    
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
            # Open FastA file    
            if filename.endswith(('.fasta', '.fna', '.fa', '.fasta.gz', '.fna.gz', '.fa.gz')):
                yield from SeqIO.parse(file, 'fasta')
            else:  
                raise ValueError(f"file format not supported: '{filename}'")
            
def kmers(sequence, k)-> dict[str:tuple[int:int]]:
        '''
        Collects all possible on the sequence passed in parameter associated with their occurrences.

        Parameters:
        :param query: a nucleotide sequence corresponding to the query in which we are looking for every Kmers possible

        Returns:
        :return dict[str:tuple[int:int]]: A dictionary with Kmers sequences as keys and a list of occurrences positions
        on the query as values.

        Examples:
        >>> a = Kmers(3)
        >>> a.query_kmers('ATCGAAATCG')
        {'ATC': [(0, 2), (6, 8)], 'TCG': [(1, 3), (7, 9)], 'CGA': [(2, 4)], 'GAA': [(3, 5)], 'AAA': [(4, 6)], 'AAT': [(5, 7)]}
        
        # Cas extrême : séquence trop courte
        >>> a.query_kmers('A')  
        {}
        '''
        #Initializes a default dictionary with list as values
        kmers = defaultdict(list)

        #Loop through every position that can contain a kmer on the query sequence
        for pos in range(len(sequence)-k+1):
            #Extract the kmers and add it as a key to the dictionary and then add their occurrences positions as values 
            kmers[sequence[pos:pos+k]].append((pos,pos+k-1))

        #Transform the default dictionary into a Python dictionary
        return dict(kmers)
    
if __name__ == '__main__':
    k = 'CCC'
    for s in read_gz("Level0.fa.gz"):
        seq = str(s.seq)
        bf = bloom_filter(seq, 3)
        if k in bf:
            print("oui")
        else:
            print("c mor")
