#Import des modules nécessaires
from collections import defaultdict
from typing import Dict, List, Set, Tuple

class DeBruijnGraph:
    """
    Classe représentant un graphe de de Bruijn basé sur un dictionnaire de k-mers.
    Elle permet d'extraire les contigs de manière sûre et efficace.
    """
    ## Intialisation du graphe 
    def __init__(self, kmers_dict: Dict[str, int]):
        self.kmers_dict = kmers_dict
        self.__graph = defaultdict(list)
        self.__reverse_graph = defaultdict(list)
        self.__build_graph()

    def __build_graph(self) -> None:
        """Construit les connexions entre préfixes et suffixes de k-mers."""
        for kmer in self.kmers_dict:
            if len(kmer) < 2:
                continue
            prefix = kmer[:-1]
            suffix = kmer[1:]
            self.__graph[prefix].append(suffix)
            self.__reverse_graph[suffix].append(prefix)

    def reconstruct_graph(self):
        self.__graph.clear()
        self.__reverse_graph.clear()
        self.__build_graph()

    ## Méthode get
    def get_successors(self, node: str) -> List[str]:
        """Retourne la liste des successeurs du nœud donné."""
        return self.__graph.get(node, [])

    def get_predecessors(self, node: str) -> List[str]:
        """Retourne la liste des prédécesseurs du nœud donné."""
        return self.__reverse_graph.get(node, [])
    
    def get_graph(self):
        return self.__graph
    
    def get_reverse_graph(self):
        return self.__reverse_graph
    
    ## Construction de contigs
    def __simple_path(self, start_node: str) -> Tuple[List[str], Set[str]]:
        """
        Construit un chemin simple à partir d'un nœud de départ et retourne les k-mers utilisés.
        """
        path = [start_node]
        used_kmers = set()

        # Extension vers l'avant
        current = start_node
        while True:
            successors = self.get_successors(current)
            if len(successors) != 1:
                break
            next_node = successors[0]
            kmer = current + next_node[-1]
            if kmer not in self.kmers_dict:
                break
            path.append(next_node)
            used_kmers.add(kmer)
            current = next_node

        # Extension vers l'arrière
        current = start_node
        while True:
            predecessors = self.get_predecessors(current)
            if len(predecessors) != 1:
                break
            prev_node = predecessors[0]
            kmer = prev_node + current[-1]
            if kmer not in self.kmers_dict:
                break
            path.insert(0, prev_node)
            used_kmers.add(kmer)
            current = prev_node

        return path, used_kmers
    
    def __assemble_sequence(self, path: List[str]) -> str:
        """Assemble une séquence à partir d'un chemin de nœuds."""
        if not path:
            return ""
        return path[0] + ''.join(node[-1] for node in path[1:])
    
    def get_all_contigs(self, output_file: str = "output_file.fa", tip_threshold = 3, bubble_threshold = 50) -> None:
        """
        Extrait tous les contigs du graphe, les écrit dans un fichier FASTA, 
        et supprime les k-mers utilisés.
        """
        self.remove_tips(tip_threshold)
        self.remove_bubbles(bubble_threshold)

        contig_num = 1
        
        with open(output_file, 'w') as out:
            for start_kmer in list(self.kmers_dict.keys()):
                if start_kmer not in self.kmers_dict:
                    continue  # K-mer déjà traité

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
    
    ## Option d'assemblage
    # Gestion des tips
    def is_tip(self, path, threshold):
        if not path:
            return False

        return len(self.get_successors(path[-1])) == 0 and len(path) < threshold
    
    def find_all_tips(self, threshold =10) -> List[List[str]]:
        tips = []

        visited = set()

        for node in self.__graph:
            if node in visited:
                continue

            path, used_kmers = self.__simple_path(node)

            visited.update(path)

            if self.is_tip(path, threshold):
                tips.append((path, used_kmers))

        return tips

    def remove_tips(self, threshold: int = 2):
        """
        Remove tips from the graph

        Parameter:
        threshold: Maximum allowed tip length

        Examples:

        # Cas de base : une branche principale avec un court tip
        #>>> g = DeBruijnGraph({'ATG':1, 'TGG':1, 'GGA':1, 'TGT':1})
        #>>> g.remove_tips(2)
        #>>> dict(g.get_graph())
        #{'AT': ['TG'], 'TG': ['GG'], 'GG': ['GA']}

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
    
        for tip in tips:
            path = self.__simple_path(tip)

            # On vérifie que le chemin se termine, sans successeur
            if len(path) <= threshold and len(self.get_successors(path[-1])) == 0:

                # Supprimer tous les k-mers correspondant à ce chemin
                for i in range(len(path)-1):
                    kmer = path[i] + path[i+1][-1]
                    if kmer in self.kmers_dict:
                        del self.kmer_dict[kmer]

        self.reconstruct_graph()

    # Gestion des bulles
    def find_all_bubbles(self, max_length= 50):
        bubbles = []
        visited = set()

        for node in self.__graph:
            if node in visited:
                continue

            successors = self.get_successors(node)
            if len(successors) < 2:
                continue  # Pas de bulle possible

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

    def remove_bubbles(self, threshold: int =50):
        """
        Remove bubbles from the graph

        Parameter:
        threshold: Maximum path length to consider when detecting bubbles

        Examples:

        # Cas de base : deux chemins partant d'un même point et se rejoignant
        >>> g = DeBruijnGraph({'ATG':1, 'TGG':1, 'TGA':1, 'GGT':1, 'GAT':1})
        >>> g.remove_bubbles()
        >>> sorted(dict(g.get_graph()).items())
        [('AT', ['TG']), ('GG', ['GT']), ('TG', ['GG'])]

        # Cas limite : les deux chemins font la même longueur
        >>> g = DeBruijnGraph({'ATG':1, 'TGA':1, 'TGC':1, 'GAA':1, 'GCA':1})
        >>> g.remove_bubbles()
        >>> sorted(dict(g.get_graph()).items())
        [('AT', ['TG']), ('TG', ['GC']), ('GC', ['A'])]

        # Cas limite : les chemins divergents ne convergent pas
        >>> g = DeBruijnGraph({'ATG':1, 'TGA':1, 'TGC':1})
        >>> g.remove_bubbles()
        >>> sorted(dict(g.get_graph()).items())
        [('AT', ['TG']), ('TG', ['GA', 'GC'])]

        # Cas complexe : deux bulles se rejoignent en un point non commun
        >>> g = DeBruijnGraph({'ATG':1, 'TGA':1, 'TGC':1, 'GAC':1, 'GCT':1})
        >>> g.remove_bubbles()
        >>> sorted(dict(g.get_graph()).items())
        [('AT', ['TG']), ('GA', ['C']), ('GC', ['T']), ('TG', ['GA', 'GC'])]
        """
        # Liste de tuples(noeuds, chemins)
        bubbles = self.find_all_bubbles(threshold) 

        for node, paths in bubbles:
            if len(paths) <2:

                # Pas de bulle si un seul chemin
                continue

            # Tous les autres sont supprimés
            for path in paths[1:]:
                if len(path) <= threshold:
                    for i in range(len(path)-1):
                        kmer = path[i] + path[i + 1][-1]
                        if kmer in self.kmers_dict:
                            del self.kmers_dict[kmer]

        self.reconstruct_graph()

