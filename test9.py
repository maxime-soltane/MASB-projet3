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

    def simple_path(self, kmer: str) -> List[str]:
        path = []

        start = kmer[:-1]
        if start not in self.graph:
            return path

    # Backward: chercher le début du chemin
        current = start
        while True:
            preds = self.get_predecessors(current)
            if len(preds) != 1:
                break
            if len(self.get_successors(preds[0])) != 1:
                break
            current = preds[0]
    
        path.append(current)

        while True:
            succs = self.get_successors(current)
            if len(succs) != 1:
                break
            if len(self.get_predecessors(succs[0])) != 1:
                break
            current = succs[0]
            path.append(current)

        return path


    def assemble_sequence(self, kmer) -> str:
        """
        Assembles a sequence from the simple path.

        Returns:
        An assembled nucleotidic sequence

        Example:
        >>> graph = DeBruijnGraph({'ATG': 1, 'TGG': 1, 'GGT': 1})
        >>> graph.assemble_sequence()
        'ATGGT'
        """
        path = self.simple_path(kmer)

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
    assembled_seq = dbg.assemble_sequence("CCCACGGACGCCAGAACGGGC")
    print("\nLong sequence test:")
    print(f"Original length: {len(longer_sequence)}")
    print(f"Assembled length: {len(assembled_seq)}")
    print(f"Correct assembly: {longer_sequence == assembled_seq}")
