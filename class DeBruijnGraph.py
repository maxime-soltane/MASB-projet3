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
        self.__ant = ant if ant is not None else []
        self.__post = post if post is not None else []

    def get_seq(self):
        return self.__sequence
    
    def get_ant(self):
        return self.__ant
    
    def get_post(self):
        return self.__post

    def add_post(self, post):
        self.__post.append(post)
    
    def add_ant(self, ant):
        self.__ant.append(ant)

    def __str__(self):
        
        # Récupère les séquences des noeuds antécédents
        ant_seqs = [node.get_seq() for node in self.get_ant()] if self.get_ant() else []
    
        # Récupère les séquences des noeuds successeurs
        post_seqs = [node.get_seq() for node in self.get_post()] if self.get_post() else []
    
        return (f"Node(sequence='{self.get_seq()}', "
                f"antecedents={ant_seqs}, "
                f"successeurs={post_seqs})")
    
class DeBruijnGraph : 
    
    def __init__(self, kmers_dict):
        self.__nodes = {} #self.__create_nodes()  
        self.__kmers_dict = kmers_dict      
        self.__build_graph()
        self.__first_node = self.__get_first_node()

    def __get_first_node (self):
        most_common_kmer = max(self.__kmers_dict, key = self.__kmers_dict.get)
        occurences = self.__kmers_dict[most_common_kmer]
        print (f"the most common node is {most_common_kmer} with {occurences} occurrences : it's the selected first node.")
        return self.__nodes[most_common_kmer]

    def __create_all_nodes(self):
        for kmer in self.__kmers_dict:
            if kmer not in self.__nodes:
                self.__nodes[kmer] = Node(kmer)

    def __nodes_connections (self):
        for kmer, count in self.__kmers_dict.items():
            suff = kmer[1:]
            current_node = self.__nodes[kmer]

            for nuc in "ATCG":
                next_kmer = suff + nuc
                if next_kmer in self.__nodes:
                    next_node = self.__nodes[next_kmer]
                    for _ in range(count):
                        current_node.add_post(next_node)
                        next_node.add_ant(current_node)

    def __build_graph(self):
        self.__create_all_nodes()
        self.__nodes_connections()

    def simple_path(self, kmer):
        if kmer not in self.__nodes:
            return []
        
        longest_path = []
        #Stack the elements : current node, current path and visited node <- to avoid cycle
        path_to_explore = [(self.__nodes[kmer], [kmer], {kmer})]
        
        
        while path_to_explore:
            current_node, path, visited = path_to_explore.pop()

            if len(path) > len(longest_path):
                    longest_path = path
                
            post = current_node.get_post()
            for node in post:
                node_seq = node.get_seq()
                if node_seq not in visited:
                    new_visited = visited.copy()
                    new_visited.add(node_seq)
                    new_path = path + [node_seq]
                    path_to_explore.append((node, new_path, new_visited))

        return longest_path

    def assemble_sequence(self):
        path = self.simple_path()
        if not path:
            return ""
        sequence = path[0]
        for node in path[1:]:
            sequence += node[-1]  # on ajoute seulement la dernière lettre
        return sequence
    
if __name__ == "__main__":
    f = read_gz("Level0.fa.gz")
    kmers_dict = {}
    for seq in f:
        kmers_dict.update(kmers(str(seq.seq), 5))
    g = DeBruijnGraph(kmers_dict)
    print(g.simple_path("CCCAC"))

#CCCACGGACGCCAGAACGGGCGTTCTCCCTAGCGTGCGCCCTGCAGAACGTTCGCGAGAACGACAGAACTCACGGACGTTCTCCCTATCGACCGTGCGCAAGAACGTCCGGCCGTACGCCCTATAGAACGAGCGCCCGCTCGGCCGTGCTATAGAACTCTCGGCCTCACGGAAGAACGTTCGTGCGCCAGAACTATCTCACGCCCTAAAGTG