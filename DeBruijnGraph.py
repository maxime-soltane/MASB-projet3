from collections import defaultdict
from typing import Dict, List, Set, Tuple

class Graph:
    """
    Classe représentant un graphe de de Bruijn basé sur un dictionnaire de k-mers.
    Elle permet d'extraire les contigs de manière sûre et efficace.
    """

    def __init__(self, kmers_dict: Dict[str, int]):
        self.kmers_dict = kmers_dict
        self.graph = defaultdict(list)
        self.reverse_graph = defaultdict(list)
        self.__build_graph()

    def __build_graph(self) -> None:
        """Construit les connexions entre préfixes et suffixes de k-mers."""
        for kmer in self.kmers_dict:
            if len(kmer) < 2:
                continue
            prefix = kmer[:-1]
            suffix = kmer[1:]
            self.graph[prefix].append(suffix)
            self.reverse_graph[suffix].append(prefix)

    def get_successors(self, node: str) -> List[str]:
        """Retourne la liste des successeurs du nœud donné."""
        return self.graph.get(node, [])

    def get_predecessors(self, node: str) -> List[str]:
        """Retourne la liste des prédécesseurs du nœud donné."""
        return self.reverse_graph.get(node, [])

    def __get_path_and_kmers(self, start_node: str) -> Tuple[List[str], Set[str]]:
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

    def get_all_contigs(self, output_file: str) -> None:
        """
        Extrait tous les contigs du graphe, les écrit dans un fichier FASTA, 
        et supprime les k-mers utilisés.
        """
        contig_num = 1
        with open(output_file, 'w') as out:
            for start_kmer in list(self.kmers_dict.keys()):
                if start_kmer not in self.kmers_dict:
                    continue  # K-mer déjà traité

                start_node = start_kmer[:-1]
                path, used_kmers = self.__get_path_and_kmers(start_node)

                if not path:
                    continue

                for k in used_kmers:
                    self.kmers_dict.pop(k, None)

                contig = self.__assemble_sequence(path)
                out.write(f">contig_{contig_num}_len_{len(contig)}\n")
                out.write('\n'.join(contig[i:i+80] for i in range(0, len(contig), 80)) + '\n')
                contig_num += 1

        print(f"Contigs générés : {contig_num - 1}, Kmers restants : {len(self.kmers_dict)}")


    def is_tip (self, path, threshold = 1):
        if not path:
            return False

        return len(self.get_successors(path[-1])) == 0 and len(path) < threshold
    
    def is_bubble(self, path, max_length = 50):
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
        
    
    
    def _find_all_bubbles(self, max_length):
        pass

    def complex_pattern_removing (self):
        pass
