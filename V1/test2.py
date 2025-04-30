from collections import defaultdict
from typing import Dict, List, Set, Tuple

class Graph:
    def __init__(self, kmers_dict: Dict[str, int]):
        self.kmers_dict = kmers_dict  # On conserve la référence originale
        self.graph = defaultdict(list)
        self.reverse_graph = defaultdict(list)
        self._build_graph()
    
    def _build_graph(self) -> None:
        """Construit le graphe une seule fois"""
        for kmer in self.kmers_dict:
            if len(kmer) < 2:
                continue
            prefix, suffix = kmer[:-1], kmer[1:]
            self.graph[prefix].append(suffix)
            self.reverse_graph[suffix].append(prefix)
    
    def get_all_contigs(self, output_file: str) -> None:
        """Version optimisée en mémoire qui supprime les kmers traités"""
        contig_num = 1
        
        with open(output_file, 'w') as out:
            # On travaille sur une copie des clés pour éviter les modifications pendant l'itération
            for start_kmer in list(self.kmers_dict.keys()):
                if start_kmer not in self.kmers_dict:  # Déjà traité
                    continue
                    
                start_node = start_kmer[:-1]
                path, used_kmers = self._get_path_and_kmers(start_node)
                
                if not path:
                    continue
                    
                # Suppression explicite des kmers utilisés
                for k in used_kmers:
                    self.kmers_dict.pop(k, None)
                
                # Assemblage du contig
                contig = path[0] + ''.join([node[-1] for node in path[1:]])
                out.write(f">contig_{contig_num}_len_{len(contig)}\n")
                out.write('\n'.join([contig[i:i+80] for i in range(0, len(contig), 80)]) + '\n')
                contig_num += 1
    
    def _get_path_and_kmers(self, start_node: str) -> Tuple[List[str], Set[str]]:
        """Retourne le chemin ET les kmers utilisés"""
        path = [start_node]
        used_kmers = set()
        
        # Extension vers l'avant
        current = start_node
        while True:
            next_nodes = self.graph.get(current, [])
            if len(next_nodes) != 1:
                break
            next_node = next_nodes[0]
            kmer = current + next_node[-1]
            if kmer not in self.kmers_dict:  # Kmer déjà traité
                break
            path.append(next_node)
            used_kmers.add(kmer)
            current = next_node
        
        # Extension vers l'arrière
        current = start_node
        while True:
            prev_nodes = self.reverse_graph.get(current, [])
            if len(prev_nodes) != 1:
                break
            prev_node = prev_nodes[0]
            kmer = prev_node + current[-1]
            if kmer not in self.kmers_dict:  # Kmer déjà traité
                break
            path.insert(0, prev_node)
            used_kmers.add(kmer)
            current = prev_node
            
        return path, used_kmers