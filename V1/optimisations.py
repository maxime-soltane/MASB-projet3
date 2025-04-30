import gzip
from Bio import SeqIO   
from DeBruijnGraph import * 
from collections import defaultdict
from time import time
from typing import Dict

def read_gz(filename: str):
        """
        Open a Fasta or FastQ file, compressed (.gz) or not, in reading mode and return sequence's iterator.

        Parameters:
        filename: the name of the file to work on

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
            
def kmers (sequence: str, k: int) -> Dict[str, int]:
    """
    Generates all kmers from a sequence and count their occurrences.

    Parameters:
    sequence: A nucleotidic sequence
    k: The size of the kmer

    Returns:
    A dictionary of kmers and their counts

    Examples:
    >>> kmers("ATCGGCAT", 3)
    {'ATC': 1, 'TCG': 1, 'CGG': 1, 'GGC': 1, 'GCA': 1, 'CAT': 1}
    >>> kmers("AAAAAA", 2)
    {'AA': 5}
    >>> kmers("ATG", 5)
    {}
    """
    kmers = defaultdict(int)
    #parcours toutes els positions où des kmers sont possibles
    for pos in range(len(sequence)-k+1):
        #les ajoute au dico avec leurs nombre d'occurrences
        kmers[sequence[pos:pos+k]] += 1
    return kmers

def kmers_filter(kmers_dict: Dict[str, int], threshold: int= 1) -> Dict[str, int]:
    """
    Filters kmers according to a minimum count threshold.

    Parameters:
    kmers_dict: The dictionary of kmer counts
    threshold: The minimum count to keep a kmer

    Returns:
    A filtered dictionnary of kmers

    Examples:
    >>> k_dict = {"A" : 4, "B" : 1, "C" : 12, "D" : 5}
    >>> kmers_filter(k_dict, 5)
    {'C': 12, 'D': 5}
    >>> kmers_filter({'A': 2, 'T': 3}, 3)
    {'T': 3}
    """
    #initialise un dico de kmer filtrer avec leur occurrences
    f_kmers = defaultdict(int)
    for kmer, count in kmers_dict.items():
        if count >= threshold:
            f_kmers[kmer] = count
    
    return f_kmers

class DeBruijnGraph:
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
        current = kmer[:-1]
    
        if current not in self.graph:
            return [], set()

        while True:

            preds = self.get_predecessors(current)
            if len(preds) != 1 or len(self.get_successors(preds[0])) != 1:
                break
                #possiblement un bulle ici donc ajouter bubble_removing ici

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