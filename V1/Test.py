from collections import defaultdict, deque
from typing import Dict, List, Set, Tuple

class Graph:
    def __init__(self, kmers_dict: Dict[str, int]):
        self.kmers_dict = kmers_dict.copy()
        self.graph = defaultdict(list)
        self.reverse_graph = defaultdict(list)
        self._build_graph()
    
    def _build_graph(self) -> None:
        """Safe graph construction with kmer validation"""
        for kmer in self.kmers_dict:
            if len(kmer) < 2:
                continue
            prefix, suffix = kmer[:-1], kmer[1:]
            self.graph[prefix].append(suffix)
            self.reverse_graph[suffix].append(prefix)
    
    def assemble_contigs(self, output_file: str = "contigs.fasta", min_length: int = 100) -> int:
        """
        Robust assembly with duplicate prevention and length filtering
        Returns number of contigs written
        """
        contig_count = 0
        used_nodes = set()
        
        with open(output_file, 'w') as out:
            # Process nodes in consistent order
            for start_node in sorted(self.graph.keys()):
                if start_node in used_nodes:
                    continue
                    
                path = self._extend_path(start_node)
                
                if not path or len(path) < 2:
                    continue
                    
                # Check if entire path is unused
                if any(node in used_nodes for node in path):
                    continue
                    
                used_nodes.update(path)
                contig = self._path_to_sequence(path)
                
                if len(contig) >= min_length:
                    out.write(f">contig_{contig_count+1}_len_{len(contig)}\n")
                    out.write('\n'.join([contig[i:i+60] for i in range(0, len(contig), 60)]) + '\n')
                    contig_count += 1
                    
        return contig_count
    
    def _extend_path(self, start_node: str) -> List[str]:
        """Extends path maximally in both directions"""
        path = [start_node]
        
        # Extend forward
        current = start_node
        while True:
            next_nodes = self.graph.get(current, [])
            if len(next_nodes) != 1:
                break
            next_node = next_nodes[0]
            if len(self.reverse_graph.get(next_node, [])) != 1:
                break
            path.append(next_node)
            current = next_node
        
        # Extend backward
        current = start_node
        while True:
            prev_nodes = self.reverse_graph.get(current, [])
            if len(prev_nodes) != 1:
                break
            prev_node = prev_nodes[0]
            if len(self.graph.get(prev_node, [])) != 1:
                break
            path.insert(0, prev_node)
            current = prev_node
            
        return path
    
    def _path_to_sequence(self, path: List[str]) -> str:
        """Converts path to DNA sequence"""
        if not path:
            return ""
        return path[0] + ''.join(node[-1] for node in path[1:])