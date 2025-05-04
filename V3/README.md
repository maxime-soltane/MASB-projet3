# MASB-projet3 : Assembleur par graphe de De Bruijn

Modules nécessaires au bon fonctionnement du script :
- Biopython (pour 'SeqIO')
- matplotlib (optionnel, pour la visualisation des histogrammes d'abondance de kmers)
- argparse
- gzip
- collections (defaultdict)

L'objectif de ce projet est de coder un assembleur basé sur la structure du graphe de De Bruijn, afin de reconstruire des contigs à partir de reads. Différents datasets présentant des particularités spécifiques devront pouvoir être assemblés pour attester de la robustesse de l'assembleur. Pour en évaluer la qualité, on soumettra les fichiers de sortie à l'analyse de l'outil QUAST. Ce projet comprend la lecture de fichiers Fasta, Fastq, qu'ils soient compressés (.gz) ou non. De plus, il permet l'extraction des kmers, la construction du graphe, la détection et gestion de structures parasites (tips et bulles). Enfin, il générera un fichier Fasta contenant les contigs assemblés.

Afin d'y parvenir, une classe à été produite :
- DBG (DeBruijn Graph)

Cette dernière, ses méthodes ainsi que le script python sont fonctionnels.

# Script

Le fichier script permet l'ouverture et la lecture de fichier aux formats Fasta, Fastq, avec ou sans compression .gz

- read_gz(filename): Retourne un itérateur sur les séquences contenues dans un fichier
- kmers(sequence, k): Extrait tous les kmers de taille k d'une séquence
- kmers_filter(kmer_dict, threshold): Filtre les kmers avec une abondance inférieure au threshold
- abundance_hist(kmers_dict): Génère un histogramme log de la distribution des kmers

# Classe DBG

La classe DBG modélise le graphe de De Bruijn, gère les tips ainsi que les bulles.

Constructeur:

- __init__(kmers_dict): Initialise le graphe à partir d'un dictionnaire de kmers

Méthodes principales:

- get_graph(): Retourne le graphe direct
- get_reverse_graph(): Retourne le graphe inverse
- get_successors(node): Retourne les successeurs d'un noeud
- get_predecessors(node): Retourne les prédécesseurs d'un noeud
- get_all_contigs(): Reconstruit les contigs et les écrit dans un fichier Fasta
- __assemble_sequence(path): Assemble une séquence à partir d'un chemin
- __simple_path(start_node): Extrait un chemin simple depuis un noeud donné

Gestion des tips:

- is_tip(path, threshold): Renvoie True si le chemin est un tip
- find_all_tips(threshold): Retourne la liste de tous les tips détectés
- remove_tips(threshold): Supprime les tips du graphe et des kmers

Gestion des bulles:

- find_all_bubbles(max_length): détecte les bulles divergentes-convergentes
- remove_bubbles(threshold): Supprime les bulles les plus courtes

# Programme principal (Main.py)

L'exécution principale se fait via le fichier "Main.py". Il permet d'exécuter l'ensemble du processus en ligne de commande. Il lit les fichiers, extrait et filtre les kmers, affiche un histogramme d'abondance si demandé, construit le graphe et génère le fichier des contigs. Le module argparse est utilisé afin d'ajouter une liste d'arguments au programme principal.

Le fichier de sortie contient un ou plusieurs contigs au format Fasta. Pour chaque contig, deux éléments sont fournis:

- Identifiant du contig: sous la forme >contig_numéro_len_longueur, où:
    - numéro est un identifiant unique pour chaque contig
    - longueur est la longueur de la sequence assemblée
- Séquence du contig: La séquence nucléotidique reconstituée à partir du graphe

Liste des arguments:

Arguments obligatoires:

- -r, --reads_file: Fichier contenant les reads
- -k, --kmers_length: Taille des kmers à extraire

Arguments optionnels:

- -o, --outfile: Fichier Fasta de sortie contenant les contigs (défaut: "output_file.fa")
- -kf, --kmers_filter_threshold: Seuil minimal d'abondance des kmers à conserver
- -tt, --tip_threshold: seuil maximal pour considérer un chemin comme un tip
- -a, --assembler: Lance l'assemblage
- -kh, --kmers_abundance_hist: Génère l'histogramme d'abondance des kmers

Exemples d'utilisation:

    - python3 Main.py -r Level1.fa.gz -o contigs_level1.fa -k 21 -a
Produira le fichier "contigs_Level1.fa" contenant les contigs assemblés à partir des reads du fichier "Level1.fa.gz" en utilisant une taille de kmers de 21.

    - python3 Main.py -r Level2.fa.gz -o contigs_Level2.fa -k 21 -kf 2 -a
Produira le fichier "contigs_Level2.fa" contenant les contigs assemblés à partir des reads du fichier "Level2.fa.gz" en utilisant une taille de kmers de 21 et en filtrant les kmers pour ne garder que ceux dont l'abondance est supérieure à 2.

    - python3 Main.py -r Level3.fa.gz -o contigs_Level3.fa -k 27 -tt 5 -a
Produira le fichier "contigs_Level3.fa" contenant les contigs assemblés à partir des reads du fichier "Level3.fa.gz" en utilisant une taille de kmers de 27 et un seuil de suppression de tips de 5.

Les arguments peuvent ensuite être combinés pour configurer l'assemblage selon les critères, par exemple:

    - python3 Main.py -r Level5.fa.gz -o contigs_Level5.fa -k 31 -kf 3 -tt 4 -a
Produira le fichier "contigs_Level5.fa.gz" contenant les contigs assemblés à partir des reads du fichier "Level5.fa.gz" en utilisant une taille de kmer de 31, un seuil d'abondance de 3 pour les kmers, en filtrant les kmers pour ne garder que ceux dont l'abondance est supérieure à 3 et en utilisant un seuil de suppression des tips de 5.

# Evaluation de la qualité d'assemblage

Les contigs sont produits au format Fasta. L’évaluation peut être réalisée avec QUAST:

Exemple d'utilisation:

    - python3 quast.py -r Ref.fa contigs_Level5.fa -m 50

Liste des arguments:

- -r: Fichier Fasta contenant le génome de référence. Il servira à comparer les contigs obtenus et à calculer des statistiques de qualité telles que le taux de couverture ou le N50
- -m: Spécifie la taille minimale de la longueur des contigs (en bp) à prendre en compte.
- contig_Level5.fa: Fichier de contigs produit par l'assembleur à analyser






