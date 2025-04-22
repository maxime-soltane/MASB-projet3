from collections import defaultdict

class DeBruijnGraph :

    def __init__ (self, kmers_dict, start_node = None):
        self.__kmers_dict = kmers_dict
        self.__graph = defaultdict(list)
        self.noeud_entrant = defaultdict(int)
        self.noeud_sortant = defaultdict(int)
        self.__build_connections()
        if start_node:
            self.__start_node = start_node
        else : 
            self.__start_node = self.__find_start_node()
        

    def __build_connections (self):
        for kmer, count in self.__kmers_dict.items():
            prefixe = kmer[:-1]
            suffixe = kmer[1:]
            for _ in range(count):
                self.__graph[prefixe].append(suffixe)
                self.noeud_entrant[suffixe] +=1
                self.noeud_sortant[prefixe] += 1

    def __find_start_node(self):
        # Choisir un nœud de départ : avec un out_degree > in_degree
        for node in self.__graph:
            if self.noeud_sortant[node] > self.noeud_entrant.get(node, 0):
                return node
        # Sinon, retourner n’importe quel nœud avec des successeurs
        return next(iter(self.get_graph()))

    def get_graph (self):
        return self.__graph
    
    def get_start (self):
        return self.__start_node
    
    def simple_path(self):
        #initialisation d'une liste vide pour stocker le chemin final
        path = []

        #initialise le chemin en cours avec le noeud de départ si présent sinon une liste vide
        chemin_en_cours = [self.__start_node] if self.__start_node else []
        
        #tant qu'on a des noeuds a explorer dans le chemin temporaire
        while chemin_en_cours:
            #on définit le noeud actuel comme le dernier de chemin en cours
            current = chemin_en_cours[-1]

            #on récupère les noeuds successeurs de current dans le graph 
            successors = self.__graph[current]
            
            #Si des successeurs existent
            if successors:
                #on prend le premier successeur
                next_node = successors[0]
                #on l'ajoute au chemin en cours
                chemin_en_cours.append(next_node)
                #on le retire du graph pour éviter de le resélectionner
                self.__graph[current].remove(next_node)
            #Si on a pas de successeur
            else:
                #on a explorer tous les chemins possibles de ce noeud : on le retire du chemin en cours et on l'ajoute au chemin
                path.append(chemin_en_cours.pop())
                
        return path
    
    def assemble_sequence(self):
        path = self.simple_path()
        if not path:
            return ""

        sequence = path[-1]
        for node in reversed(path[:-1]):
            sequence += node[-1]
        
        return sequence