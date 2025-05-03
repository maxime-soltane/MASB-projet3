from collections import defaultdict
from typing import List, Set, Tuple


class DBG:

    def __init__(self, kmers_dict):
        self.__kmers_dict = kmers_dict
        self.__graph = defaultdict(list)
        self.__build_graph()
    
    def __build_graph(self):
        self.__graph = defaultdict(list)  # Réinitialise le graphe
    
        for kmer in self.__kmers_dict:
            if len(kmer) < 2:
                continue
            prefix = kmer[:-1]
            suffix = kmer[1:]
            self.__graph[prefix].append(suffix)

    #Get method

    def get_successors(self, node: str) -> List[str]:
        """
        Returns the list of successors for a given node.

        Parameter:
        node: A node in the graph

        Returns:
        A list of successors nodes
        
        Examples:
        >>> g = DBG({'ATG':1, 'TGG':1, 'GGA':1, 'TGT':1})
        >>> g.get_successors("AT")
        ['TG']
        >>> g.get_successors("TG")
        ['GG', 'GT']
        """
        return self.__graph.get(node, [])
    
    def get_predecessors (self, node: str) -> List[str]:
        """
        Returns the list of predecessors for a given node.

        Parameter:
        node: A node in the graph

        Returns:
        A list of predecessors nodes
        
        Examples:
        >>> g = DBG({'ATG':1, 'TGG':1, 'GGA':1, 'TGT':1})
        >>> g.get_predecessors("GG")
        ['TG']
        >>> g.get_predecessors("AT")
        []
        """
        predecessors = []
        for prefix, suffixes in self.__graph.items():
            if node in suffixes:
                predecessors.append(prefix)
        return predecessors

    def get_graph(self):
        """
        Returns the forward de Bruijn graph as a dictionary.

        Returns:
        A dictionary where keys are prefixes and values are lists of suffixes

        Example:
        >>> g = DBG({'ATG':1, 'TGG':1, 'GGA':1, 'TGT':1})
        >>> g.get_graph()
        {'AT': ['TG'], 'TG': ['GG', 'GT'], 'GG': ['GA']}
        """
        return dict(self.__graph)
    
    def get_kmers_dict(self):
        return self.__kmers_dict
    
    # Construction de chemin

    def __extend_forward(self, start_node):
        """
        Étend le chemin vers l'avant à partir de `start_node` tant que les nœuds ont un seul successeur.
        
        Paramètre:
        start_node: Nœud de départ.
        
        Retourne:
        Liste des nœuds parcourus (y compris `start_node`).

        >>> g = DBG({'ATG':1, 'TGG':1, 'GGA':1, 'TGT':1})
        >>> g.get_graph()
        {'AT': ['TG'], 'TG': ['GG', 'GT'], 'GG': ['GA']}
        >>> g._DBG__extend_forward('AT')
        ['AT', 'TG']
        >>> g._DBG__extend_forward('TG')
        ['TG']
        >>> g._DBG__extend_forward('GG')
        ['GG', 'GA']
        >>> g._DBG__extend_forward('GT')
        ['GT']
        >>> g._DBG__extend_forward('GG')
        ['GG', 'GA']
        >>> g._DBG__extend_forward('GA')
        ['GA']
        """
        path = [start_node]
        current_node = start_node

        while True:
            successors = self.get_successors(current_node)
            if len(successors) != 1:
                break  # Bifurcation ou fin du chemin

            next_node = successors[0]
            predecessors_of_next = self.get_predecessors(next_node)

            if len(predecessors_of_next) != 1 or predecessors_of_next[0] != current_node:
                break  # Successeur a plusieurs prédécesseurs ou lien non-unique

            path.append(next_node)
            current_node = next_node

        return path

    def __extend_backward(self, start_node):
        """
        Étend le chemin vers l'arrière à partir de `start_node` tant que les nœuds ont un seul prédécesseur.
        
        Paramètre:
        start_node: Nœud de départ.
        
        Retourne:
        Liste des nœuds parcourus (sans `start_node`, pour éviter les duplications).

        >>> g = DBG({'ATG':1, 'TGG':1, 'GGA':1, 'TGT':1})
        >>> g.get_graph()
        {'AT': ['TG'], 'TG': ['GG', 'GT'], 'GG': ['GA']}
        >>> g._DBG__extend_backward('AT')
        []
        >>> g._DBG__extend_backward('TG')
        ['AT']
        >>> g._DBG__extend_backward('GG')
        []
        >>> g._DBG__extend_backward('GT')
        []
        >>> g._DBG__extend_backward('GG')
        []
        >>> g._DBG__extend_backward('GA')
        ['GG']
        """
        path = []
        current_node = start_node

        while True:
            predecessors = self.get_predecessors(current_node)
            if len(predecessors) != 1:
                break  # Bifurcation ou début du chemin

            pred = predecessors[0]
            successors_of_pred = self.get_successors(pred)

            if len(successors_of_pred) != 1 or successors_of_pred[0] != current_node:
                break  # Prédécesseur a plusieurs successeurs ou lien non-unique

            path.insert(0, pred)  # On ajoute au début pour garder l'ordre
            current_node = pred

        return path

    def __simple_path(self, start_node):
        """
        Combine l'extension vers l'arrière et l'avant pour former un chemin linéaire maximal.
        
        Paramètre:
        start_node: Nœud de départ.
        
        Retourne:
        Chelin linéaire complet (arrière + avant).

        >>> g = DBG({'ATG':1, 'TGG':1, 'GGA':1, 'TGT':1})
        >>> g.get_graph()
        {'AT': ['TG'], 'TG': ['GG', 'GT'], 'GG': ['GA']}
        >>> g._DBG__simple_path('AT')
        ['AT', 'TG']
        >>> g._DBG__simple_path('TG')
        ['AT', 'TG']
        >>> g._DBG__simple_path('GG')
        ['GG', 'GA']
        >>> g._DBG__simple_path('GT')
        ['GT']
        >>> g._DBG__simple_path('GA')
        ['GG', 'GA']
        """
        backward_path = self.__extend_backward(start_node)
        forward_path = self.__extend_forward(start_node)
        full_path = backward_path + forward_path  # Fusion sans duplication de `start_node`

        return full_path


    # Options d'assemblage

    def is_tip(self, path: List[str], threshold: int=5) -> bool:
        """
        Checks if a given path is a tip (a short dead-end path).

        Parameters:
        path: The list of nodes representing a path
        threshold: The maximum length to be considered a tip

        Returns:
        True if it's a tip, False otherwise

        Examples:
        >>> g = DBG({'ATG':1, 'TGG':1, 'GGC':1})
        >>> p = g._DBG__simple_path('AT')
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
        >>> g = DBG({'ATG':1, 'TGG':1, 'GGC':1})
        >>> g._DBG__simple_path('AT')
        ['AT', 'TG', 'GG', 'GC']
        >>> g.find_all_tips()
        [['AT', 'TG', 'GG', 'GC']]

        #Cas 2 : 3 tips
        >>> g = DBG({'ATG': 1, 'TGG': 1, 'GGA': 1, 'GGT': 1, 'GTC': 1, 'TCC': 1, 'TCA': 1})
        >>> g.get_graph()
        {'AT': ['TG'], 'TG': ['GG'], 'GG': ['GA', 'GT'], 'GT': ['TC'], 'TC': ['CC', 'CA']}
        >>> g.find_all_tips(2)
        [['GG', 'GA'], ['TC', 'CC'], ['TC', 'CA']]

        >>> # Cas 3: Aucun tip (tous les chemins ont des successeurs)
        >>> g = DBG({'ATG': 1, 'TGG': 1, 'GGG': 1, 'GGT': 1, 'GTT': 1, 'TTT': 1})
        >>> g.get_graph()
        {'AT': ['TG'], 'TG': ['GG'], 'GG': ['GG', 'GT'], 'GT': ['TT'], 'TT': ['TT']}
        >>> len(g.find_all_tips())
        0

        # Cas 4: Tip plus long que le seuil (ne doit pas être détecté)
        >>> g = DBG({'ATG':1, 'TGC':1, 'GCA':1, 'CAA':1})
        >>> g.get_graph()
        {'AT': ['TG'], 'TG': ['GC'], 'GC': ['CA'], 'CA': ['AA']}
        >>> len(g.find_all_tips(3))  
        0
        """
        tips = []
        visited = set()

        for node in self.__graph:
            if node in visited:
                continue

            path = self.__simple_path(node)
            if self.is_tip(path, threshold):
                tips.append(path)
                visited.update(path)
            
            successors = self.get_successors(node)
            if len(successors) > 1:
                for succ in successors:
                    if succ in visited:
                        continue
                    path = self.__simple_path(succ)
                    if self.is_tip(path, threshold):
                        full_path = [node] + path
                        tips.append(full_path)
                        visited.update(path)

        print("tips trouvés")
        return tips

    def remove_tips(self, threshold=5):
        """
        Remove tips from the graph based on a length threshold.
 
        Parameter:
        threshold: The maximum length to be considered a tip
 
        Parameter:
        threshold: Maximum allowed tip length
 
        Examples:
        # Cas de base : une branche principale avec un court tip
        >>> g = DBG({'ATG':1, 'TGG':1, 'GGA':1, 'TGT':1})
        >>> g.get_graph()
        {'AT': ['TG'], 'TG': ['GG', 'GT'], 'GG': ['GA']}
        >>> g.remove_tips(2)
        >>> g.get_graph()
        {'AT': ['TG'], 'TG': ['GG'], 'GG': ['GA']}
        >>> g.get_kmers_dict()
        {'ATG': 1, 'TGG': 1, 'GGA': 1}
        """
        k = len(next(iter(self.__kmers_dict))) 
        tips = self.find_all_tips(threshold)
        
        for tip in tips:  
                
            # Reconstruit la séquence et extrait les k-mers
            sequence = tip[0] + ''.join(node[-1] for node in tip[1:])
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                self.__kmers_dict.pop(kmer, None)
        
        self.__build_graph()
        print("tips éliminés")

# Bubbles management

    def remove_bubbles(self) -> List[Tuple[str, List[List[str]]]]:
        """
        Detects all bubbles in the graph and remove paths to let 1 path reamining.

        Returns:
        A list of tuples containing a starting node and a list of alternative paths that converge later

        >>> kmers_bulle = {'ATG':1, 'TGC':1, 'GCA':1, 'CAA':1, 'GCT':1, 'CTA':1, 'TAA':1}
        >>> g = DBG(kmers_bulle)
        >>> g.get_graph()
        {'AT': ['TG'], 'TG': ['GC'], 'GC': ['CA', 'CT'], 'CA': ['AA'], 'CT': ['TA'], 'TA': ['AA']}
        >>> g.remove_bubbles()
        >>> g.get_graph()
        {'AT': ['TG'], 'TG': ['GC'], 'GC': ['CA'], 'CA': ['AA']}
        """
        k = len(next(iter(self.__kmers_dict))) 
        
        #pour chaque noeud
        for node in self.__graph:
            #on collecte ses successeurs
            successors = self.get_successors(node)
            #si plus d'un successeur
            if len(successors) >= 2:
                #point de convergence possible = 
                convergence_point = []
                #pour chaque successeur dans successeur
                for su in successors:
                    #on construit le chemin du successseur actuel
                    path = self.__simple_path(su)
                    #on prend le dernier noeud du chemin
                    last_node = path[-1]
                    #on récupère ses successeurs
                    s = self.get_successors(last_node)
                    #s'il en a qu'un = point de convergence
                    if len(s) == 1:
                        #si le successeur n'est pas ajouté dans point de convergence possible
                        if s[0] not in convergence_point:
                            #on l'ajoute
                            convergence_point.append(s[0])
                        #si dans convergence point = point de convergence de la bulle
                        elif s[0] in convergence_point:
                            #on ajoute au chemin le noeud de départ de la bulle
                            path.insert(0, node)
                            #on ajoute le point de convergence
                            path.append(s[0])
                            #on assemble la séquence présente dans le chemin
                            sequence = self.__assemble_sequence(path)
                            #on reconstruit les kmers et les enlève du dico
                            for i in range(len(sequence) - k + 1):
                                kmer = sequence[i:i+k]
                                self.__kmers_dict.pop(kmer, None)

        self.__build_graph()
        print("bulles éliminées")

# Assemblage de séquence

    def __assemble_sequence(self, path: List[str]) -> str:
        """
        Assembles a sequence from a nodes-made path.
        
        Parameter:
        path: A list of nodes representing the contig path
        
        returns:
        A single string representing the assembled sequence

        Examples:
        #Cas de base
        >>> g = DBG({'ATG':1, 'TGG':1, 'GGC':1})        
        >>> g._DBG__assemble_sequence(['AT', 'TG', 'GG', 'GC'])
        'ATGGC'

        #Cas 2 : chemin vide
        >>> g = DBG({})        
        >>> g._DBG__assemble_sequence([])
        ''
        """
        if not path:
            return ""
        return path[0] + ''.join(node[-1] for node in path[1:])

    def get_all_contigs(self, output_file: str = "output_file.fa", tip_threshold = 3) -> None:
        """
        Extracts all the contigs from the graph, writes them in a fasta file and deletes the used kmers.

        Parameters:
        output_file: Path to the output fasta file
        tip_threshold: The threshold for tip removal
        bubble_threshold: The threshold for bubble removal   
        """
        self.remove_tips(tip_threshold)
        self.remove_bubbles()

        contig_num = 1
        k = len(next(iter(self.__kmers_dict))) 

        with open(output_file, 'w') as out:
            for start_kmer in list(self.__kmers_dict.keys()):
                if start_kmer not in self.__kmers_dict:
                    # Kmers already processed
                    continue  

                start_node = start_kmer[:-1]
                path = self.__simple_path(start_node)

                if not path:
                    continue

                sequence = self.__assemble_sequence(path)
                for i in range(len(sequence) - k + 1):
                                kmer = sequence[i:i+k]
                                self.__kmers_dict.pop(kmer, None)

                contig = self.__assemble_sequence(path)
                out.write(f">contig_{contig_num}_len_{len(contig)}\n")
                out.write('\n'.join(contig[i:i+60] for i in range(0, len(contig), 60)) + '\n')
                contig_num += 1

        print(f"Contigs générés : {contig_num - 1}")
        print(f"{output_file} was generated\n")