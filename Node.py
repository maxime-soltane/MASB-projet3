class Node:

# on commence du kmer le plus vu mais pas dans la vrai vie -> prenons la médiane 
#vérifier si on a des N
    def __init__(self, sequence, ant: list =  None, post: list = None):
        self.__sequence = sequence
        self.__ant = ant
        self.__post = post

class Graphe_de_de_Bruijn:

    def __init__(self, kmers_dict):
        self.__nodes_list = []
        self.__kmers_dict = kmers_dict

    def find_nodes(self, node: Node):
        suff = node.__sequence[-1]
        if self.__kmers_dict.keys().startwith(suff+"A"):
            n = Node(kmer_seq_on_ne_sait_pas_comment, ant = node)
            node.__post.append(n)
        elif self.__kmers_dict.keys().startwith(suff+"T"):
            n = Node(kmer_seq_on_ne_sait_pas_comment, ant = node)
            node.__post.append(n)
        elif self.__kmers_dict.keys().startwith(suff+"C"):
            n = Node(kmer_seq_on_ne_sait_pas_comment, ant = node)
            node.__post.append(n)
        elif self.__kmers_dict.keys().startwith(suff+"G"):
            n = Node(kmer_seq_on_ne_sait_pas_comment, ant = node)
            node.__post.append(n)
        self.__nodes_list.append()