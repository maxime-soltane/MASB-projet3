from collections import defaultdict
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

class DeBruijnGraph :

    def __init__ (self, kmers_dict, start_node = None):
        self.__kmers_dict = kmers_dict
        self.__graph = defaultdict(list)
        self.noeud_entrant = defaultdict(int)
        self.noeud_sortant = defaultdict(int)
        self.__build_connections()
        if start_node:
            self.__start_node = start_node
        else : 
            self.__start_node = self.__find_start_node()
        

    def __build_connections (self):
        for kmer, count in self.__kmers_dict.items():
            prefixe = kmer[:-1]
            suffixe = kmer[1:]
            for _ in range(count):
                self.__graph[prefixe].append(suffixe)
                self.noeud_entrant[suffixe] +=1
                self.noeud_sortant[prefixe] += 1

    def __find_start_node(self):
        # Choisir un nœud de départ : avec un out_degree > in_degree
        for node in self.__graph:
            if self.noeud_sortant[node] > self.noeud_entrant.get(node, 0):
                return node
        # Sinon, retourner n’importe quel nœud avec des successeurs
        return next(iter(self.get_graph()))

    def get_graph (self):
        return self.__graph
    
    def get_start (self):
        return self.__start_node
    
    def simple_path(self):
        path = []

        stack = [self.__start_node] if self.__start_node else []
        
        while stack:
            current = stack[-1]
            successors = self.__graph.get(current, [])
            
            if successors:
                next_node = successors[0]
                stack.append(next_node)
                self.__graph[current].remove(next_node)
            else:
                path.append(stack.pop())
        return path
    
    def assemble_sequence(self):
        path = self.simple_path()
        if not path:
            return ""

        sequence = path[-1]
        for node in reversed(path[:-1]):
            sequence += node[-1]
        
        return sequence
    
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