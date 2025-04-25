from collections import defaultdict
from typing import Dict, List
class Graph:
    """
    A class that represents a de Bruijn graph built from a collection of kmers.
    This graph will be used to assemble genomes.
    """

    def __init__(self, kmers_dict: dict):
        self.__kmers_dict = kmers_dict
        self.graph = defaultdict(list)
        self.reverse_graph = defaultdict(list)
        self.__build_connections()

    def __build_connections (self) -> None:
        """
        Builds graph connections from the kmers, each kmer is split into a prefixe and a suffix.
        """
        for kmer in self.__kmers_dict:
            prefix = kmer[:-1]
            suffix = kmer[1:]
            self.graph[prefix].append(suffix)
            self.reverse_graph[suffix].append(prefix)

    def get_graph(self) -> Dict[str, List[str]]:
        """
        Gets the graph representation.

        Returns:
        A dictionary of prefix
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
        >>> graph = Graph({'ATG':1})
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
        >>> graph = Graph({'ATG':1})
        >>> graph.get_predecessors('TG')
        ['AT']
        """
        return self.reverse_graph.get(node, [])
    
    def __simple_path(self, kmer: str) -> List[str]:
        """
        Computes a simple path starting from a given kmer.

        Parameters:
        kmer: A starting kmer from the kmer dictionary

        Returns:
        A list of nodes constructing the path

        Example:
        >>> g = Graph({'ATG':1, 'TGG':1, 'GGA':1})
        >>> g._Graph__simple_path('TGG')
        ['AT', 'TG', 'GG', 'GA']
        """
        path = []
        start = kmer[:-1]
    
        if start not in self.graph:
            return path

        current = start  
        while True:
            preds = self.get_predecessors(current)
            if len(preds) != 1:
                break
            if len(self.get_successors(preds[0])) != 1:
                break

            kmer_to_remove = preds[0] + current[-1]
            if kmer_to_remove in self.__kmers_dict:
                del self.__kmers_dict[kmer_to_remove]
        
            current = preds[0]
    
        path.append(current)

        while True:
            succs = self.get_successors(current)
            if len(succs) != 1:
                break
            if len(self.get_predecessors(succs[0])) != 1:
                break
        
            kmer_to_remove = current + succs[0][-1]
            if kmer_to_remove in self.__kmers_dict:
                del self.__kmers_dict[kmer_to_remove]

            current = succs[0]
        
            path.append(current)

        return path
    
    def __assemble_sequence(self, kmer: int) -> str:
        """
        Assembles a contig sequence starting from the given kmer.

        Parameter:
        kmer: A starting kmer

        Returns:
        A string representing the assembled contig

        Example:
        >>> graph = Graph({'ATG':1, 'TGG':1, 'GGA':1})
        >>> graph._Graph__assemble_sequence('TGG')
        'ATGGA'
        """
        path = self.__simple_path(kmer)

        if not path:
            return ""
        
        contig = path[0]
        
        for node in path[1:]:
            contig += node[-1]

        return contig
    

    def get_all_contigs(self, output_file: str = "output.fasta") -> None:
        """
        Assemble all contigs and write them to a Fasta file.

        Parameter:
        output_file: output Fasta filename
        """
        contig_count = 1
        with open(output_file, 'w') as f:
            while self.__kmers_dict:
                
                kmer, _ = self.__kmers_dict.popitem()

                contig = self.__assemble_sequence(kmer)

                if contig:
                    f.write(f">contig_{contig_count}_from_{kmer}\n")
                    for i in range(0, len(contig), 60):
                        f.write(contig[i:i+60] + '\n')
                    contig_count += 1

        print("Le fichier a été créé.")