#Import des modules nécessaires
from collections import defaultdict
from typing import Dict, List, Set, Tuple

class DeBruijnGraph:
    """
    A Class that represents a de Bruijn graph made from a dictionary of kmers.
    It allows safe and effcicient extraction of contigs.

    """
    ## Graph's initialization
    def __init__(self, kmers_dict: Dict[str, int]) -> None:
        self.kmers_dict = kmers_dict
        self.__graph = defaultdict(list)
        self.__reverse_graph = defaultdict(list)
        self.__build_graph()

    def __build_graph(self) -> None:
        """
        Builds the graph by linking kmer prefixes and suffixes.
        """
        for kmer in self.kmers_dict:
            if len(kmer) < 2:
                continue
            prefix = kmer[:-1]
            suffix = kmer[1:]
            self.__graph[prefix].append(suffix)
            self.__reverse_graph[suffix].append(prefix)

    def reconstruct_graph(self) -> None:
        """
        Clears and rebuilds the graph by using the current kmer dictionary.
        """
        self.__graph.clear()
        self.__reverse_graph.clear()
        self.__build_graph()

    ## Get method
    def get_successors(self, node: str) -> List[str]:
        """
        Returns the list of successors for a given node.

        Parameter:
        node: A node in the graph

        Returns:
        A list of successors nodes
        
        Examples:
        >>> g = DeBruijnGraph({'ATG':1, 'TGG':1, 'GGA':1, 'TGT':1})
        >>> g.get_successors("AT")
        ['TG']
        >>> g.get_successors("TG")
        ['GG', 'GT']
        """
        return self.__graph.get(node, [])

    def get_predecessors(self, node: str) -> List[str]:
        """
        Returns the list of predecessors for a given node.

        Parameter:
        node: A node in the graph

        Returns:
        A list of predecessors nodes
        
        Examples:
        >>> g = DeBruijnGraph({'ATG':1, 'TGG':1, 'GGA':1, 'TGT':1})
        >>> g.get_predecessors("GG")
        ['TG']
        >>> g.get_predecessors("AT")
        []
        """
        return self.__reverse_graph.get(node, [])
    
    def get_graph(self) -> Dict[str, List[str]]:
        """
        Returns the forward de Bruijn graph as a dictionary.

        Returns:
        A dictionary where keys are prefixes and values are lists of suffixes

        Example:
        >>> g = DeBruijnGraph({'ATG':1, 'TGG':1, 'GGA':1, 'TGT':1})
        >>> g.get_graph()
        {'AT': ['TG'], 'TG': ['GG', 'GT'], 'GG': ['GA']}
        """
        return dict(self.__graph)
    
    def get_reverse_graph(self) -> Dict[str, List[str]]:
        """
        Returns the reverse de Bruijn graph.

        Returns:
        A dictionary where keys are suffixes and values are lists of prefixes

        Example:
        >>> g = DeBruijnGraph({'ATG':1, 'TGG':1, 'GGA':1})
        >>> g.get_reverse_graph()
        {'TG': ['AT'], 'GG': ['TG'], 'GA': ['GG']}
        """
        return dict(self.__reverse_graph)
    
    ## Construction de contigs
    def __simple_path(self, start_node: str) -> Tuple[List[str], Set[str]]:
        """
        Builds a simple path from a starting node and returns the used kmers.

        Parameter:
        start_node: The starting node of the path

        Returns:
        A tuple containing the path as a list of nodes and the set of kmers used

        Examples:
        #Cas de base 
        >>> g = DeBruijnGraph({'ATG':1, 'TGG':1, 'GGC':1})
        >>> g.get_graph()
        {'AT': ['TG'], 'TG': ['GG'], 'GG': ['GC']}
        >>> g._DeBruijnGraph__simple_path('AT')[0]
        ['AT', 'TG', 'GG', 'GC']

        #Cas 2 : Gestion d'un cycle
        >>> g = DeBruijnGraph({'ATG':1, 'TGG':1, 'GGA':1, 'GAT':1})
        >>> g.get_graph()
        {'AT': ['TG'], 'TG': ['GG'], 'GG': ['GA'], 'GA': ['AT']}
        >>> g._DeBruijnGraph__simple_path('AT')[0]
        ['AT', 'TG', 'GG', 'GA']

        #Cas 3 : Start_node absent dans le graphe
        >>> g = DeBruijnGraph({'ATG':1, 'TGG':1, 'GGA':1, 'GAT':1})
        >>> g.get_graph()
        {'AT': ['TG'], 'TG': ['GG'], 'GG': ['GA'], 'GA': ['AT']}
        >>> g._DeBruijnGraph__simple_path('TT')
        ([], set())
        """
        if not start_node in self.get_graph():
            return [], set()
        
        path = [start_node]
        used_kmers = set()

        # Forward extension
        current = start_node
        while True:
            successors = self.get_successors(current)
            if len(successors) != 1:
                break
            next_node = successors[0]
            if next_node in path:
                break
            kmer = current + next_node[-1]
            if kmer not in self.kmers_dict:
                break
            path.append(next_node)
            used_kmers.add(kmer)
            current = next_node

        # Backward extension
        current = start_node
        while True:
            predecessors = self.get_predecessors(current)
            if len(predecessors) != 1:
                break
            prev_node = predecessors[0]
            if prev_node in path:
                break 
            kmer = prev_node + current[-1]
            if kmer not in self.kmers_dict:
                break
            path.insert(0, prev_node)
            used_kmers.add(kmer)
            current = prev_node

        return path, used_kmers
    
    def __assemble_sequence(self, path: List[str]) -> str:
        """
        Assembles a sequence from a nodes-made path.
        
        Parameter:
        path: A list of nodes representing the contig path
        
        returns:
        A single string representing the assembled sequence

        Examples:
        #Cas de base
        >>> g = DeBruijnGraph({'ATG':1, 'TGG':1, 'GGC':1})        
        >>> g._DeBruijnGraph__assemble_sequence(['AT', 'TG', 'GG', 'GC'])
        'ATGGC'

        #Cas 2 : chemin vide
        >>> g = DeBruijnGraph({})        
        >>> g._DeBruijnGraph__assemble_sequence([])
        ''
        """
        if not path:
            return ""
        return path[0] + ''.join(node[-1] for node in path[1:])
    
    def get_all_contigs(self, output_file: str = "output_file.fa", tip_threshold = 3, bubble_threshold = 50) -> None:
        """
        Extracts all the contigs from the graph, writes them in a fasta file and deletes the used kmers.

        Parameters:
        output_file: Path to the output fasta file
        tip_threshold: The threshold for tip removal
        bubble_threshold: The threshold for bubble removal   
        """
        self.remove_tips(tip_threshold)
        self.remove_bubbles(bubble_threshold)

        contig_num = 1
        
        with open(output_file, 'w') as out:
            for start_kmer in list(self.kmers_dict.keys()):
                if start_kmer not in self.kmers_dict:
                    # Kmers already processed
                    continue  

                start_node = start_kmer[:-1]
                path, used_kmers = self.__simple_path(start_node)

                if not path:
                    continue

                for k in used_kmers:
                    self.kmers_dict.pop(k, None)

                contig = self.__assemble_sequence(path)
                out.write(f">contig_{contig_num}_len_{len(contig)}\n")
                out.write('\n'.join(contig[i:i+60] for i in range(0, len(contig), 60)) + '\n')
                contig_num += 1

        print(f"Contigs générés : {contig_num - 1}")
        print(f"{output_file} was generated\n")
    
    # Asembly options

    # Tips management

    def is_tip(self, path: List[str], threshold: int=5) -> bool:
        """
        Checks if a given path is a tip (a short dead-end path).

        Parameters:
        path: The list of nodes representing a path
        threshold: The maximum length to be considered a tip

        Returns:
        True if it's a tip, False otherwise

        Examples:
        >>> g = DeBruijnGraph({'ATG':1, 'TGG':1, 'GGC':1})
        >>> p = g._DeBruijnGraph__simple_path('AT')[0]
        >>> print(p)
        ['AT', 'TG', 'GG', 'GC']
        >>> g.is_tip(p)
        True
        >>> g.is_tip(["AT", "TG", "GG", "GC"], 2)
        False
        """
        if not path or len(path) >= threshold:
            return False
        
        return len(self.get_successors(path[-1])) == 0
    
    def find_all_tips(self, threshold:int=5) -> List[Tuple[List[str], Set[str]]]:
        """
        Detects all tips in the graph under a certain threshold

        Parameter:
        threshold: The maximum path length considered as a tip

        Returns:
        A list of tuples containing the path and the kmers used

        Examples:
        #Cas 1 : 1 seul tip
        >>> g = DeBruijnGraph({'ATG':1, 'TGG':1, 'GGC':1})
        >>> p = g._DeBruijnGraph__simple_path('AT')[0]
        >>> g.find_all_tips()[0][0] 
        ['AT', 'TG', 'GG', 'GC']

        #Cas 2 : 2 tips
        >>> g = DeBruijnGraph({'ATG':1, 'TGG':1, 'GGC':1, 'CAA':1, 'AAT':1, 'GTT':1, 'TTG':1})
        >>> print(g.get_graph())
        {'AT': ['TG'], 'TG': ['GG'], 'GG': ['GC'], 'CA': ['AA'], 'AA': ['AT'], 'GT': ['TT'], 'TT': ['TG']}
        >>> tips = g.find_all_tips(10)
        >>> len(tips)
        2

        >>> # Cas 3: Aucun tip (tous les chemins ont des successeurs)
        >>> g = DeBruijnGraph({'ATG':1, 'TGT':1, 'GTG':1, 'TGG':1, 'GGT':1, 'GTG':1})
        >>> g.get_graph()
        {'AT': ['TG'], 'TG': ['GT', 'GG'], 'GT': ['TG'], 'GG': ['GT']}
        >>> len(g.find_all_tips())
        0

        >>> # Cas 4: Tip plus long que le seuil (ne doit pas être détecté)
        >>> g = DeBruijnGraph({'ATG':1, 'TGG':1, 'GGC':1, 'GCC':1, 'CCA':1})
        >>> path, _ = g._DeBruijnGraph__simple_path('AT')
        >>> print(path)
        ['AT', 'TG', 'GG', 'GC', 'CC', 'CA']
        >>> len(g.find_all_tips(threshold=3))  
        0
        """
        tips = []
        visited = set()

        for node in sorted(self.__graph.keys()):
            if node in visited:
                continue

            path, used_kmers = self.__simple_path(node)

            if self.is_tip(path, threshold):
                tips.append((path, used_kmers))
                visited.update(path)
            elif len(path) > 1:
                visited.update(path[:-1])
            else:
                visited.add(node)

        return tips
        
    def remove_tips(self, threshold: int = 5) -> None:
        """
        Remove tips from the graph based on a length threshold.

        Parameter:
        threshold: The maximum length to be considered a tip

        Parameter:
        threshold: Maximum allowed tip length

        Examples:
        # Cas de base : une branche principale avec un court tip
        >>> g = DeBruijnGraph({'ATG':1, 'TGG':1, 'GGA':1, 'TGT':1})
        >>> g.get_graph()
        {'AT': ['TG'], 'TG': ['GG', 'GT'], 'GG': ['GA']}
        >>> g.remove_tips(2)
        >>> g.get_graph()
        {'AT': ['TG'], 'TG': ['GG', 'GT'], 'GG': ['GA']}

        # Cas limite : faux positif (le tip a un successeur)
        #>>> g = DeBruijnGraph({'ATG':1, 'TGG':1, 'GGA':1, 'TGT':1, 'GTT':1})
        #>>> g.remove_tips(2)
        #>>> dict(g.get_graph())
        #{'AT': ['TG'], 'TG': ['GG', 'GT'], 'GG': ['GA'], 'GT': ['TT']}

        # Cas complexe : deux tips symétriques
        #>>> g = DeBruijnGraph({'ATG':1, 'TGG':1, 'GGG':1, 'TGA':1, 'GAC':1})
        #>>> g.remove_tips(2)
        #>>> sorted(dict(g.get_graph()).items())
        #[('AT', ['TG']), ('GG', ['GG']), ('TG', ['GG'])]
        """
        tips = self.find_all_tips(threshold)
    
        for path, used_kmers in tips:

            # Check if the path ends, without successors
            if len(path) <= threshold and len(self.get_successors(path[-1])) == 0:

                # Delete all kmers corresponding to this path
                for kmer in used_kmers:
                    if kmer in self.kmers_dict:
                        del self.kmers_dict[kmer]

        self.reconstruct_graph()

# Bubbles management

    def find_all_bubbles(self, max_length= 50) -> List[Tuple[str, List[List[str]]]]:
        """
        Detects all bubbles in the graph under a maximum path length.

        Parameter:
        max_length: The maximum path length to explore when searching for bubbles

        Returns:
        A list of tuples containing a starting node and a list of alternative paths that converge later
        """
        bubbles = []
        visited = set()

        for node in self.__graph:
            if node in visited:
                continue

            successors = self.get_successors(node)
            if len(successors) < 2:
                # No bubble possible
                continue 

            convergence_points = defaultdict(list)

            for succ in successors:
                path = [node, succ]
                current = succ

                while len(path) <= max_length:
                    next_nodes = self.get_successors(current)
                    if len(next_nodes) != 1:
                        break

                    current = next_nodes[0]
                    path.append(current)

                    if len(self.get_predecessors(current)) > 1 and self.get_successors(current) == 1:
                        convergence_points[current].append(path.copy())
                        break

            for conv_node, paths in convergence_points.items():
                if len(paths) >= 2:
                    bubbles.append((node, paths))
                    visited.update(p for path in paths for p in path)

        return bubbles

    def remove_bubbles(self, threshold: int =50) -> None:
        """
        Remove bubbles from the graph based on a length threshold.

        Parameter:
        threshold: Maximum path length to consider when detecting bubbles

        Examples:

        # Cas de base : deux chemins partant d'un même point et se rejoignant
        >>> g = DeBruijnGraph({'ATG':1, 'TGG':1, 'TGA':1, 'GGT':1, 'GAT':1})
        >>> g.remove_bubbles()
        >>> sorted(dict(g.get_graph()).items())
        [('AT', ['TG']), ('GA', ['AT']), ('GG', ['GT']), ('TG', ['GG', 'GA'])]

        # Cas limite : les deux chemins font la même longueur
        >>> g = DeBruijnGraph({'ATG':1, 'TGA':1, 'TGC':1, 'GAA':1, 'GCA':1})
        >>> g.remove_bubbles()
        >>> sorted(dict(g.get_graph()).items())
        [('AT', ['TG']), ('GA', ['AA']), ('GC', ['CA']), ('TG', ['GA', 'GC'])]

        # Cas limite : les chemins divergents ne convergent pas
        >>> g = DeBruijnGraph({'ATG':1, 'TGA':1, 'TGC':1})
        >>> g.remove_bubbles()
        >>> sorted(dict(g.get_graph()).items())
        [('AT', ['TG']), ('TG', ['GA', 'GC'])]

        # Cas complexe : deux bulles se rejoignent en un point non commun
        >>> g = DeBruijnGraph({'ATG':1, 'TGA':1, 'TGC':1, 'GAC':1, 'GCT':1})
        >>> g.remove_bubbles()
        >>> sorted(dict(g.get_graph()).items())
        [('AT', ['TG']), ('GA', ['AC']), ('GC', ['CT']), ('TG', ['GA', 'GC'])]
        """
        # List of tuples(nodes, paths)
        bubbles = self.find_all_bubbles(threshold) 

        for node, paths in bubbles:
            if len(paths) <2:

                # No bubble if only one path
                continue

            # Every others are deleted
            for path in paths[1:]:
                if len(path) <= threshold:
                    for i in range(len(path)-1):
                        kmer = path[i] + path[i + 1][-1]
                        if kmer in self.kmers_dict:
                            del self.kmers_dict[kmer]

        self.reconstruct_graph()
        
if __name__ == '__main__':
    g = DeBruijnGraph({'ATG':1, 'TGG':1, 'GGC':1, 
                       'CAA':1, 'AAT':1, 'GTT':1, 
                       'TTG':1, 'TGA':1})
    print(g.get_graph())
    print(g._DeBruijnGraph__simple_path('AT'))
    print(g._DeBruijnGraph__simple_path('GT'))
    print(g._DeBruijnGraph__simple_path('CA'))

    

