class Node:

    # on commence du kmer le plus vu mais pas dans la vrai vie -> prenons la médiane 
    #vérifier si on a des N
    def __init__(self, sequence, ant: list =  None, post: list = None):
        self.__sequence = sequence
        self.__ant = ant
        self.__post = post

        