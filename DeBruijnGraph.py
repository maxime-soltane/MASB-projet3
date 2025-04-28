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
        (['AT', 'TG', 'GG', 'GA'], {'TGG', 'ATG', 'GGA'})
        """
        path = []
        kmers_in_path = set()
        start = kmer[:-1]
    
        if start not in self.graph:
            return path, kmers_in_path

        current = start  
        while True:
            preds = self.get_predecessors(current)
            if len(preds) != 1 or len(self.get_successors(preds[0])) != 1:
                break

            kmer_to_remove = preds[0] + current[-1]
            kmers_in_path.add(kmer_to_remove)

            current = preds[0]
    
        path.append(current)

        while True:
            succs = self.get_successors(current)
            if len(succs) != 1 or len(self.get_predecessors(succs[0])) != 1:
                break
        
            kmer_to_remove = current + succs[0][-1]
            kmers_in_path.add(kmer_to_remove)

            current = succs[0]
        
            path.append(current)

        return path, kmers_in_path
    
    def __assemble_sequence(self, path) -> str:
        """
        Assembles a contig sequence starting from the given kmer.

        Parameter:
        kmer: A starting kmer

        Returns:
        A string representing the assembled contig

        Example:
        >>> graph = Graph({'ATG':1, 'TGG':1, 'GGA':1})
        >>> path, _ = graph._Graph__simple_path('TGG')
        >>> graph._Graph__assemble_sequence(path)
        'ATGGA'
        """
        if not path:
            return ""
        
        return path[0] + "".join(node[-1] for node in path[1:])
    
    def get_all_contigs(self, output_file: str = "output.fasta") -> None:
        """
        Assemble all contigs and write them to a Fasta file.

        Parameter:
        output_file: output Fasta filename
        """
        contig_count = 1
        with open(output_file, 'w') as f:
            while self.__kmers_dict:
                
                kmer = next(iter(self.__kmers_dict))

                path, kmers_in_path = self.__simple_path(kmer)
                
                cleaned_path, kmers_to_remove = self.tip_removal(path, kmers_in_path)

                if not cleaned_path:
                    for kmer in kmers_to_remove: 
                        self.__kmers_dict.pop(kmer, None)

                contig = self.__assemble_sequence(cleaned_path)

                for kmer in kmers_in_path:
                    self.__kmers_dict.pop(kmer, None)
                
                f.write(f">contig_{contig_count}_of_length_{len(contig)}\n")
                for i in range(0, len(contig), 60):
                    f.write(contig[i:i+60] + '\n')
                contig_count += 1

        print(f"Contigs written to {output_file}")
    
    def tip_removal(self, path, kmers_in_path, threshold=1):
        """
        >>> graph = Graph({'ATG':1, 'TGG':1, 'GGA':1})
        >>> path, kmers_in_path = graph._Graph__simple_path('TGG')
        >>> graph.tip_removal(path, kmers_in_path)
        (['AT', 'TG', 'GG', 'GA'], set())
        >>> graph.tip_removal(path, kmers_in_path, 4)
        ([], { 'ATG', 'TGG', 'GGA'})
        >>> graph.tip_removal([], {})
        ([], set())
        """
        #chemin vide
        if not path :
            return path, set()
    
        #si court et dernier sans successeur = dead-end = tip
        if not self.get_successors(path[-1])and len(path) <= threshold:
            return [], kmers_in_path
        
        return path, set()