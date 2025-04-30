#Import des modules nécessaires au main
import argparse
from collections import defaultdict
from time import time

#Import des fonctions pour la lecture de fichier et l'extraction de kmers
from Script import *

#Import des méthodes pour créer le Graphe de DeBruijn
from DeBruijnGraph import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #Paramètres des fichiers
    parser.add_argument("-r", "--reads_file", required=True, type = str, 
                        help = "Reads the assembler will work on")
    parser.add_argument("-o", "--outfile", required = False, type = str, 
                        help = "Output file name with contigs abtained with the assembler")
    
    #Paramètres liés au kmers
    parser.add_argument("-k", "--kmers_length", required = True, type = int, 
                        help = "Length of kmers to extract")
    parser.add_argument("-kf", "--kmers_filter_threshold", required=False, type = int, 
                        help = "Abundance minimal of kmers for them to being kept")
    parser.add_argument("-kh", "--kmers_abundance_hist", required= False, action='store_true', 
                        help="Construct Kmers abundance histogram")
    
    #Paramètres liés à l'assembleur
    parser.add_argument("-a", "--assembler", required=False, action='store_true',
                        help = "Si présent assemble")
    parser.add_argument("-tt", "--tip_threshold", required=False, type=int,
                        help = "Max length of an alternative path to be considered as a tip")
    parser.add_argument("-bt", "--bubble_threshold", required=False, type=int,
                        help = "Max length of an alternative path to be considered as a bubble")

    args = parser.parse_args()

    start = time()
    f = read_gz(args.reads_file)

    kmers_dict = defaultdict(int)
    for seq in f:
        for kmer, count in kmers(str(seq.seq), args.kmers_length).items():
            kmers_dict[kmer] += count

    #Création de l'histogramme d'abondance des kmers
    if args.kmers_abundance_hist:
        abundance_hist(kmers_dict)

    if args.kmers_filter_threshold:
        f_kmers = kmers_filter(kmers_dict, args.kmers_filter_threshold)
        kmers_dict = f_kmers

    print("Kmers dictionnary generated")

    if args.assembler:
        dbg = DeBruijnGraph(kmers_dict)
        print("DeBruijn graph generated")
        
        #Gestion des arguments optionnels
        #3 arguments présents
        if args.outfile and args.tip_threshold and args.bubble_threshold:
            dbg.get_all_contigs(args.outfile, args.tip_threshold, args.bubble_threshold)

        #2 arguments présents
        elif args.outfile and args.tip_threshold:
                dbg.get_all_contigs(args.outfile, args.tip_threshold)
        elif args.outfile and args.bubble_threshold:
                dbg.get_all_contigs(args.outfile, args.bubble_threshold)
        elif args.tip_threshold and args.bubble_threshold: 
            dbg.get_all_contigs(args.tip_threshold, args.bubble_threshold)
        
        #1 argument présent
        elif args.outfile: 
            dbg.get_all_contigs(args.outfile)
        elif args.tip_threshold:
            dbg.get_all_contigs( args.tip_threshold)
        elif args.bubble_threshold: 
            dbg.get_all_contigs(args.bubble_threshold)
        
        #0 argument présent
        else:
            dbg.get_all_contigs()

        end = time()
        print(f"Execution time : {end-start}\n")