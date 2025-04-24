from collections import defaultdict
from typing import Dict, List
from Script import read_gz
from Script import kmers

class DeBruijnGraph :
    """
    A class that represents a de Bruijn graph built with a collection of kmers.
    The graph is used to assembly genomes by overlapping kmers. 
    """

    def __init__(self, kmers_dict: Dict[str, int]) -> None:
        self.__kmers_dict = kmers_dict
        self.graph = defaultdict(list)
        self.__build_connections()

    def __build_connections (self) -> None:
        """
        Builds graph connections from the kmers, each kmer is split into a prefixe and a suffix.
        """
        for kmer in self.__kmers_dict:
            prefix = kmer[:-1]
            suffix = kmer[1:]
            self.graph[prefix].append(suffix)

    def get_graph(self) -> Dict[str, List[str]]:
        """
        Gets the graph representation.
        """
        return self.graph
    
    def get_successors(self, node: str) -> List[str]:
        """
        Gets the list of successors for a node.

        Parameters:
        node: A node in the graph

        Returns:
        A list of successor nodes

        Example:
        >>> graph = DeBruijnGraph({'ATG':1})
        >>> graph.get_successors('AT')
        ['TG']
        """
        return self.graph.get(node, [])

    def get_predecessors(self, node: str) -> List[str]:
        """
        Gets the list of predecessors for a node.

        Parameters:
        node: A node in the graph

        Returns:
        A list of predecessors nodes

        Example:
        >>> graph = DeBruijnGraph({'ATG':1})
        >>> graph.get_predecessors('TG')
        ['AT']
        """
        preds = []
        for prefix in self.graph:
            for suffix in self.graph[prefix]:
                if suffix == node:
                    preds.append(prefix)
        return preds

    def simple_path(self, kmer) -> List[str]:
        """
        Tries to find a simple path starting from a node with no predecessors.
        Only follows paths with a single successor.

        Returns:
        A list of nodes forming the path

        Examples:
        >>> graph = DeBruijnGraph({'ATG':1, 'TGG':1, 'GGT':1})
        >>> graph.simple_path()
        ['AT', 'TG', 'GG', 'GT']

        # No clear start node
        >>> graph = DeBruijnGraph({'AAA': 1, 'AAT': 1, 'ATA': 1, 'TAA': 1})
        >>> graph.simple_path()
        []
        """
        path = []
        #on récupère le préfixe du kmer
        pref_kmer = kmer[:-1]

        #si le préfixe du kmer n'est pas présent dans le graph le chemin est vide
        if pref_kmer not in self.get_graph().keys():
            return path
        else: #le kmer est présent: on ajoute son préfixe au chemin
            path.append(pref_kmer)

        #renvoie le suffixe
        current_node = self.get_graph()[pref_kmer] 

        #on crée une pile qui se modifie au même moment que l'avancement du graph
        pile = [current_node]

        #tant qu'on a un truc dans la pile
        while pile:
            #on l'ajoute au chemin
            path.append(pile) 
            #on récupère les successeurs du noeud actuel
            pile = self.get_successors(pile)
            #on vérifie que le noeud actuel a encore des successeurs et qu'il n'a pas plus d'un père/fils
            if len(pile) > 1 or len(pile) == 0 or self.get_predecessors(pile) > 1:
                break

        return path

    def assemble_sequence(self) -> str:
        """
        Assembles a sequence from the simple path.

        Returns:
        An assembled nucleotidic sequence

        Example:
        >>> graph = DeBruijnGraph({'ATG': 1, 'TGG': 1, 'GGT': 1})
        >>> graph.assemble_sequence()
        'ATGGT'
        """
        path = self.simple_path()

        if not path:
            return ""
        
        contig = path[0]
        
        for node in path[1:]:
            contig += node[-1]

        return contig

if __name__ == '__main__':
    # Test with the longer sequence
    f = read_gz("Level0.fa.gz")
    kmers_dict = {}
    for seq in f:
        kmers_dict.update(kmers(str(seq.seq), 21))

    longer_sequence = "CCCACGGACGCCAGAACGGGCGTTCTCCCTAGCGTGCGCCCTGCAGAACGTTCGCGAGAACGACAGAACTCACGGACGTTCTCCCTATCGACCGTGCGCAAGAACGTCCGGCCGTACGCCCTATAGAACGAGCGCCCGCTCGGCCGTGCTATAGAACTCTCGGCCTCACGGAAGAACGTTCGTGCGCCAGAACTATCTCACGCCCTAAAGTG"
    
    dbg = DeBruijnGraph(kmers_dict)
    print(dbg.simple_path("CCCACGGACGCCAGAACGGGC"))