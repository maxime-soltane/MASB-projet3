from collections import defaultdict


class DeBruijnGraph :

    def __init__(self, kmers_dict):
        self.__kmers_dict = kmers_dict
        self.graph = defaultdict(list)
        self.__build_connections()

    def __build_connections (self):
        for kmer in self.__kmers_dict:
            prefix = kmer[:-1]
            suffix = kmer[1:]
            self.graph[prefix].append(suffix)

    def get_graph(self):
        return self.graph
    
    def get_successors(self, node):
        return self.graph.get(node, [])

    def get_predecessors(self, node):
        preds = []
        for prefix in self.graph:
            for suffix in self.graph[prefix]:
                if suffix == node:
                    preds.append(prefix)
        return preds

    def simple_path(self):
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

    def assemble_sequence(self):
        path = self.simple_path()

        if not path:
            return ""
        
        contig = path[0]
        
        for node in path[1:]:
            contig += node[-1]

        return contig
