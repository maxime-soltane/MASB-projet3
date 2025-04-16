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

class Node:

    # on commence du kmer le plus vu mais pas dans la vrai vie -> prenons la médiane 
    #vérifier si on a des N
    def __init__(self, sequence, ant: list =  None, post: list = None): #ant et post sont des objets Node
        self.__sequence = sequence
        self.__ant = ant
        self.__post = post

    def get_seq(self):
        return self.__sequence
    
    def get_ant(self):
        return self.__ant
    
    def get_post(self):
        return self.__post

    def set_post(self, post):
        self.__post = post

class DeBruijnGraph : 
    
    def __init__(self, kmers_dict):
        self.__nodes_list = []  
        self.__kmers_dict = kmers_dict
        self.__first_node = self.first_node()

    def first_node (self):
        most_common_node = max(self.__kmers_dict, key = self.__kmers_dict.get)
        occurences = self.__kmers_dict[most_common_node]
        print (f"the most common node is {most_common_node} with {occurences} occurrences : it's the selected first node.")
        node = Node(most_common_node)
        self.__nodes_list.append(node)
        return node
    
    def find_nodes (self, node: Node):
        nodes_children = []
        for nuc in "ATCG":
            kmer = node.get_seq()[1:] + nuc
            if kmer in self.__kmers_dict:
                new_node = Node(kmer, ant = [node])
                self.__nodes_list.append(new_node)
                nodes_children.append(new_node)
        node.set_post(nodes_children)
        
    def build_graph(self):
        self.find_nodes(self.__first_node)
        print("Nodes in graph:")
        for node in self.__nodes_list:
            print(node)

if __name__ == "__main__":
    f = read_gz("Level0.fa.gz")
    kmer_dict = {}
    for seq in f:
        kmer_dict.update(kmers(str(seq.seq), 11))
    g = DeBruijnGraph(kmer_dict)
    g.build_graph()
