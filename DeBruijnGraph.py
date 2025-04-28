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
    
    def test(self, output_file: str = "output.fasta", tip_threshold=1) -> None:
        contig_count = 1

        with open(output_file, 'w') as f:
            while self.__kmers_dict:
                kmer = next(iter(self.__kmers_dict)) #Prend le prochain Kmers

                path, kmers_in_path = self.__simple_path(kmer)

                if self.is_tip(path, tip_threshold):
                    for kmer in kmers_in_path:
                        self.__kmers_dict.pop(kmer, None)
                    self.__build_connections()
                    continue

                else:
                    contig = self.__assemble_sequence(path)

                    for kmer in kmers_in_path:
                        self.__kmers_dict.pop(kmer, None)

                    f.write(f">contig_{contig_count}_of_length_{len(contig)}\n")
                    for i in range(0, len(contig), 60):
                        f.write(contig[i:i+60] + '\n')
                    contig_count += 1

        print(f"Contigs written to {output_file}")

    def is_tip (self, path, threshold = 1):

        if  len(self.get_successors(path[-1])) == 0:
            if len(path) < threshold:
                return True
        return False

    def get_all_contigs(self, output_file: str = "output.fasta", tip_threshold=1) -> None:
        """
        Version corrigée avec :
        - Gestion correcte des k-mers supprimés
        - Évite les doublons dans les contigs
        """
        contig_count = 1

        with open(output_file, 'w') as f:
            while self.__kmers_dict:
                kmer = next(iter(self.__kmers_dict)) #Prend le prochain Kmers

                path, kmers_in_path = self.__simple_path(kmer)

                if self.is_tip(path, tip_threshold):
                    for kmer in kmers_in_path:
                        self.__kmers_dict.pop(kmer, None)
                    self.__build_connections()
                    continue
                
                if self.is_bubble(path):
                    pass 
                    #Gestion des bulles à définir
                    
                else:
                    contig = self.__assemble_sequence(path)

                    for kmer in kmers_in_path:
                        self.__kmers_dict.pop(kmer, None)

                    f.write(f">contig_{contig_count}_of_length_{len(contig)}\n")
                    for i in range(0, len(contig), 60):
                        f.write(contig[i:i+60] + '\n')
                    contig_count += 1

        print(f"Contigs written to {output_file}")

    def is_tip (self, path, threshold = 1):

        if self.get_successors(path[-1]) == 0:
            if len(path) < threshold:
                return True
        return False
    
    def is_bubble(self, path, max_length):
        #on récupère le dernier noeud du chemin:
        last_node = path[-1]

        #on regarde s'il a plus d'un successeur
        successors = self.get_successors(last_node)

        #si moins de 2 successeurs pas une bulle
        if len(successors) <2 :
            return False
        
        #Stocker les points de convergence
        #Clé = noeud de convergence, Valeur = liste des chemins y menant
        convergence_points = defaultdict(list)

        #on fait le chemin des x chemins  = nombre de successeur
        for succ in successors:
            current_path = [last_node, succ]
            current_node = succ

            while len(current_path) <= max_length:
                next_nodes = self.get_successors(current_node)
                
                #vérifie si on a toujours qu'un chemin simple
                if len(next_nodes) != 1:
                    break
                
                current = next_nodes[0]
                current_path.append(current)

                #vérifie si ce noeud est un point de convergence
                if len(self.get_predecessors(current)) > 1 and self.get_successors(current) == 1:
                    convergence_points[current].append(path)
                    break

        #bulle existe si au moins 2 chemins convergent vers un même noeud
        return any(len(paths) >= 2 for paths in convergence_points.values())
        
    def bubble_removing (self, current_node, max_path_length = 10):
        pass


    def complex_pattern_removing (self):
        pass
