# Rapport, `Julianne Cocq, Maxime Soltane`

## __Temps d'exécution et utilisation mémoire__

### Implémentation

- Le Script :

La lecture de fichier est réalisé à l'aide de 'SeqIO' du module 'Biopython'. Celui-ci permet d'itérer sur les séquences en évitant de stocker l'ensemble des reads en mémoire.

L'utilisation de dictionnaire par défaut ("defaultdict" du module collections) dans l'ensemble du code a permis de stocker le kmers dans l'objectif de récupérer
de manière aisé et peu coûteuse les valeurs qui y sont présentes. De plus, ce type de dictionnaire permet d'éviter la vérification de présence d'une clé dans le dictionnaire.

- La classe DBG :

La construction du graphe est réalisé de telle sorte à avoir 2 graphes distincts : "graph" et "reverse_graph" qui sont des dictionnaires par défaut.
Le premier contient, comme clés, les préfixes et, en valeur associés, leurs suffixes.
Le second contient les suffixes comme clés et les préfixes qui leurs sont associés comme valeurs.
Cette implémentation de 2 graphes distincts est coûteux en mémoire (on double la mémoire utilisée) mais permet de faciliter l'accès aux prédecesseurs et successeurs des noeuds avec une complexité en O(1) (voir'get_successors' et 'get_predecessors') réduisant le coût temporelle. 

La méthode 'find_all_tips' a été optimisée en trouvant les noeuds de départ de tips potentiels avant de réaliser la vérification complète de présence de tips. En effet, les noeuds potentiels de départ de tips correspondent à soit des noeuds sans prédécesseurs (dans le cas où le tip se trouve dès le début du graphe) soit à des noeuds avec plus d'un successeur : ceci sont alors ajouter à un set (permettant d'éviter les duplications) et à partir de ce dernier, on vérifie la réelle présence de tips. Enfin, on réalise un second passage pour s'assurer que tous les noeuds ont bien été traités. 

La méthode 'get_all_contigs' écrit directement les contigs obtenus dans le fichier de sortie (évitant un stockage inutile) puis élimine les kmers déjà traités du dictionnaire pour libérer de la mémoire.  

- Conclusion :

La mise en place de structure de données comme les dictionnaires par défaut et les sets, ont permis de favoriser une gestion de la mémoire et du temps d'exécution en réduisant la complexité, facilitant l'accès et la modification des données et en évitant les opérations de vérification pour éviter les doublons. 
Enfin, certaines parties du code ont été optimisées afin d'apporter une meilleure gestion de la mémoire et du temps d'exécution. 

### Mesure du temps d'exécution et de la mémoire utilsée

Ces deux paramètres ont été mesurés à l'aide de la commande : 
```bash
"/usr/bin/time -v" python3 Main.py
```
#### Level 0

```bash
/usr/bin/time -v python3 Main.py -r Level0.fa.gz -k 21 -a
```
Maximum resident set size (kbytes): 67248

Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.78

#### Level 1

```bash
/usr/bin/time -v python3 Main.py -r Level0.fa.gz -k 21 -a
```
Maximum resident set size (kbytes): 86580

Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02.23

#### Level 2

```bash
/usr/bin/time -v python3 Main.py -r Level2.fa.gz -k 300 -kf 2 -a
```
Maximum resident set size (kbytes): 164636s

Elapsed (wall clock) time (h:mm:ss or m:ss): 0:04.28

#### Level 3

```bash
/usr/bin/time -v python3 Main.py -r Level3.fa.gz -k 800 -kf 3 -a
```
Maximum resident set size (kbytes): 969176

Elapsed (wall clock) time (h:mm:ss or m:ss): 0:11.23

#### Level 4

```bash
/usr/bin/time -v python3 Main.py -r Level4.fa.gz -k 1500 -kf 3 -tt 1000 -a
```
Maximum resident set size (kbytes): 4775344

Elapsed (wall clock) time (h:mm:ss or m:ss): 0:24.80

#### Level 5

```bash
/usr/bin/time -v python3 Main.py -r Level5.fa.gz -k 2000 -kf 3 -tt 1000 -a
```
Maximum resident set size (kbytes): 7481700

Elapsed (wall clock) time (h:mm:ss or m:ss): 1:29.93

#### Level 6

- k = 1000
```bash
/usr/bin/time -v python3 Main.py -r Level6.fa.gz -k 1000 -kf 3 -tt 1000 -a
```


- k = 2000
```bash
/usr/bin/time -v python3 Main.py -r Level5.fa.gz -k 2000 -kf 3 -tt 1000 -a
```


#### Level 7

## __Qualité d'assemblage__

La qualité d'assemblage a été évalué à l'aide de l'outil QUAST.
Les rapports obtenus se trouvent dans le dossier QUAST du projet. 

### Level 0

- Contexte

Il s'agit d'un fichier test ne comportant qu'une séquence de 212 nucléotides. 
Le QUAST a donc été utilisé uniquement pour obtenir les caractéristiques de cette séquence sans la comparer à un fichier de réference.

- Observations

On observe tout d'abord que la séquence obtenue avec l'assembleur correspond parfaitement à la séquence du fichier original de taille 212, et aucune fragmentation ni erreurs n'ont été réalisées par l'assembleur.

### Level 1

- Observations 

On observe un seul contig de 29403 nucléotides parfaitement alignés sur le génome de référence avec un 'Genome fraction' de 98,328% et un pourcentage GC similaire (37,95% contre 37,97%).
Le NGA90 vaut 29403 indique que l'unique contig couvre plus de 90% de la référence.
Aucune erreur d'assemblage et d'alignement n'est présente.
Le ratio de duplication de 1 montre que le génome de référence est bien couvert qu'une unique fois.

- Conclusion

L'assemblage est correct et sans erreur mais 500 bases sont manquantes.

### Level 2

- Observations

On obtient 2 contigs de taille 29719 et 29106 bases nucléotidiques. 
Seul le premier s'aligne sans aucune erreur sur la référence avec une couverture de 99,4%.
Le GC obtenus via assemblage vaut 39,16% contre 37,97% pour la référence : on suppose que cela est causé par le second contig non aligné. 
Le NGA90 vaut 29719 ce qui confirme que le contig aligné couvre correctement la référence.
Aucune erreur d'assemblage et d'alignement n'est présente. 
Cependant, la taille totale est de 58 825.
Le ratio de duplication de 1 montre que le génome de référence est bien couvert qu'une unique fois.

- Conclusions

Le premier contig correspond presque entièrement au génome de référence. Cependant, le second contig ne correspond pas à ce dernier.
On peut supposer qu'il s'agit d'une contamination lors du séquençage ou d'une duplicationcomme le mentionne la taille totale des contigs qui vaut environ 2 fois la taille de la référence.

### Level 3

- Observations

On a 4 contigs pour une longueur totale de 115439 bases nucléotidiques dont l'un d'entre eux est aligné avec 28833 bases. 
Le pourcentage de GC vaut 39,49% contre 37,97% pour la référence.
Le NGA90 de 28833 montre que le contig aligné couvre plus de 90% de la référence.
Le génome est couvert à hauteur de 96,422% et on a 86606 bases non alignées correspondante aux 3 autres contigs.
Le ratio de duplication de 1 montre que le génome de référence est bien couvert qu'une unique fois.
Il n'y a aucune erreur d'assemblage.

- Conclusions

Les bases non alignées représentent un problème majeur dans l'assemblage même si un contig s'aligne correctement sur le génome de référence. Ceci peut être expliquer par une gestion des erreurs incomplètes de l'assembleur.

### Level 4

- Observations 

On a 8 contigs de taille totale 230733 avec le plus long ayant 29360 bases mais 6 d'entre eux ne s'alignent pas sur la référence. 
Le taux de GC est plus élevé de 1,67%.
La référence est couverte à 96,375%.
Le NGA50 de 28789 montre que les contigs alignés couvrent plus de 50% de la référence.
Le ration de duplication est de 1,994 indique que les alignements sont couverts 2 fois.

- Conclusions

On peut déduire du ratio de duplication que les 2 contigs alignées sont similaires entraînant pratiquement un alignement doublé sur le génome de référence même si des mismatchs sont présents. 
On peut alors supposer qu'ils correspondent à de l'hétérozygotie ou une duplication lors de l'assemblage.

### Level 5

- Observations

On a 16 contigs de taille totale 455375 nucléotides dont 12 ne s'alignent pas sur la référence.
Le ratio de duplication vaut 3,870. 
Le NGA90 montre que les contigs alignés d'au moins 29154 bp couvrent 90% de la référence. En effet, il est couvert à 97,495% par les contigs alignés. 
Le taux de GC est plus élevé de 1,67%.

On a 1735 mismatches présents dans l'alignement et aucun misassemblie n'a été détecté.

- Conclusion

On a davantage de contigs non alignées mais le génome est legèrement plus couverts qu'au niveau 4. Cependant, les mismatches sont également plus nombreux. 

### Level 6

#### __k = 1000__
- Observations 

On a 94 contigs dont seulement 22 ont une taille supérieure à 25000 et dont 60 ne sont pas entiérement alignés sur la référence. 
Le taux de GC est de 39,57% dans les contigs contre 37,97% dans le génome de référence. 
Ce dernier est couvert à hauteur de 99,324% avec 7932 mismatchs. Le NGA90 montre que les contigs de taille d'au moins 28574 bp couvrent plus de 90% de la référence.
Le ratio de duplication est de 13,228.
Lors de l'alignement, on a 7932 mismatchs.

- Conclusion

On a de trop nombreux contigs expliquaient par une taille de k bien trop faible entraînant des contigs courts et créant un assemblage fragmenté. Cependant, certains contigs couvrent correctement le génome malgré des mismatchs.
Le résultat est donc un alignement de qualité faible malgré une couverture élevée.

#### __k= 2000__

- Observations

On a 32 contigs tous de taille supérieur à 25000 bp pour une taille totale de 921703 bp : 24 d'entre eux ne s'alignent pas.
Le taux de GC est légèrement plus élevé qu'avec k=1000 : 39,79%.
Le NGA90 montre que les contigs de taille d'au moins 29283 bp couvrent plus de 90 % de la référence (Genome fraction = 98,893%).
Le ratio de duplication est plus faible que précédemment et vaut 7,764.
On a 4091 mismatchs.

- Conclusion

Le taux de couverture du génome de référence est élevé mais le ratio de duplication montre de nouveau une superposition des contigs. Cependant, on peut voir qu'avec un k plus élevé, le nombre de mismatch et ce ratio diminue fortement. On a alors un alignement de meilleure qualité avec une fragmentation plus faible.

### Level 7

- Observations

On a 380 contigs dont seulement 19 de taille supérieur à 25000 bp et dont 258 complétement inalignés sur la référence. 
Le taux de GC est de 39,74%.
Le NGA90 montre que le génome de référence est couvert par les contigs de taille d'au moins 27392. Le taux de couverture vaut 98,497% mais 16183 mismatches sont présents. 
Le ratio de duplication est de 24,605.

- Conclusion

On a une fragmentation très importante malgré un taux de couverture élevé. De plus, le nombre très important de mismatches montre un problème majeur dans l'assemblage de la séquence. Ceci peut par exemple être expliqué par un k bien trop faible, cependant, l'implémentation actuelle de l'assembleur ne permet pas d'augmenter cette valeur.

## __Perspectives__

### Problèmes de gestions de Kmers

On a pu voir en expérimentant que la construction du Graphe de DeBruijn, dans son implémentation actuelle, s'exécute efficacement et rapidement. Cependant, la collecte des kmers et leur filtrage devient rapidement difficile sur les niveaux élevés rallogneant le temps d'exécution et la mémoire utilisée. On peut alors toujours construire le graphe de DeBruijn sur des Kmers de petites tailles mais y résulte un assemblage biaisé et de mauvaise qualité. Par conséquent, la fonction 'abundance_hist' ne peut plus être utilisée efficacement afin de définir le seuil d'abondance des kmers à filtrer étant donné que de trop nombreux kmers de k faible sont présents. 

### Problèmes de gestions des erreurs

Dans l'implémentation actuelle du graphe, seul les erreurs causées par les tips et les bulles sont gérées. Les répétitions et les patterns plus complexes ne sont par exemples pas traités. Ceci peut causer des assemblages biaisés avec des contigs très similaire augmentant le ratio de duplication, des msimatches et des fragmentations par exemple.

### Solutions 

L'utilisation de structure de données différentes et l'implémentation de méthodes supplémentaires pour la gestion d'erreurs d'assemblages permettraient d'obtenir de meilleurs assembleurs sur des reads plus complexes.