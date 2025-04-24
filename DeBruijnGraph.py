from collections import defaultdict
from typing import Dict, List

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

    def simple_path(self) -> List[str]:
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
        visited = set()
        for start in self.graph:
            if len(self.get_predecessors(start)) == 0:
                path = [start]
                current = start
                while current in self.graph and len(self.graph[current]) == 1:
                    next_node = self.graph[current][0]
                    if next_node in visited:
                        break
                    path.append(next_node)
                    visited.add(current)
                    current = next_node
                return path
        return []

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
