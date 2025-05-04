# Import the needed modules to the main
import argparse
from collections import defaultdict
from time import time

# Import needed functions to file reading and kmers extraction
from Script import *

# Import needed methods to create the DeBruijn graph
from DBG_V3 import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # Files parameters
    parser.add_argument("-r", "--reads_file", required=True, type = str, 
                        help = "Reads the assembler will work on")
    parser.add_argument("-o", "--outfile", required = False, type = str, 
                        help = "Output file name with contigs abtained with the assembler")
    
    # kmers parameters
    parser.add_argument("-k", "--kmers_length", required = True, type = int, 
                        help = "Length of kmers to extract")
    parser.add_argument("-kf", "--kmers_filter_threshold", required=False, type = int, 
                        help = "Abundance minimal of kmers for them to being kept")
    parser.add_argument("-kh", "--kmers_abundance_hist", required= False, action='store_true', 
                        help="Construct Kmers abundance histogram")
    
    # Assembler parameters
    parser.add_argument("-a", "--assembler", required=False, action='store_true',
                        help = "Si présent assemble")
    parser.add_argument("-tt", "--tip_threshold", required=False, type=int,
                        help = "Max length of an alternative path to be considered as a tip")

    args = parser.parse_args()

    if args.kmers_length <=0:
        raise ValueError("La taille des kmers doit être supérieure à 0.")

    if args.kmers_filter_threshold is not None and args.kmers_filter_threshold <=0:
        raise ValueError(" Le seuil pour -kf doit être strictement positif.")
    
    if args.tip_threshold is not None and args.tip_threshold <=0:
        raise ValueError(" Le seuil pour -tt doit être strictement positif.")

    start = time()
    f = read_gz(args.reads_file)

    kmers_dict = defaultdict(int)
    for seq in f:
        for kmer, count in kmers(str(seq.seq), args.kmers_length).items():
            kmers_dict[kmer] += count

    if not kmers_dict:
        raise ValueError("Aucun kmer n'a pu être extrait, vérifiez le fichier ou la valeur de k")
    
    # Creation of abundance histogram of kmers
    if args.kmers_abundance_hist:
        abundance_hist(kmers_dict)

    if args.kmers_filter_threshold:
        f_kmers = kmers_filter(kmers_dict, args.kmers_filter_threshold)
        kmers_dict = f_kmers

    print("Kmers dictionnary generated")

    if args.assembler:
        dbg = DBG(kmers_dict)
        print("DeBruijn graph generated")
        
        # Management of optional arguments
        # 2 arguments present
        if args.outfile and args.tip_threshold:
            dbg.get_all_contigs(args.outfile, args.tip_threshold)
        
        if args.outfile and not args.outfile.endswith((".fa", ".fasta")):
            raise ValueError("Le fichier doit comporter l'extension '.fa' ou '.fasta")

        # 1 argument present
        elif args.outfile: 
            dbg.get_all_contigs(args.outfile)
        elif args.tip_threshold:
            dbg.get_all_contigs(tip_threshold=args.tip_threshold)
        
        # 0 argument present
        else:
            dbg.get_all_contigs()

        end = time()
        print(f"Execution time : {end-start}\n")