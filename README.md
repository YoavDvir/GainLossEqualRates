# GainLossEqualRates

In this repository, we test methods for phylogenetic tree reconstruction for prokaryotes. For the simulations, we assume our stochastic gain/loss model.
All the methods are based on gene dynamics and equal gain and loss rates. 
The methods DCJ, GC, GCU, and CGA compute distance matrices, and then BioPython NJ reconstructs the tree.

If the DCJ method is applied, Java 8 must be installed to run. You must also download UniMoG (BiBiServ2 - DCJ), which should be in the same folder as the code.


You can find here the following 4 Python files:
1. Indel_sim.py - trees or one-edge simulations
   
This program conducts simulations:

Example for the command line: py indel_sim.py parameters_file

Example files for the parameters file that run tree simulations are:

(In these cases, it creates random trees and reconstructs them using DCJ, GC, GCU, or CGA to compute the distance matrix and then BioPython NJ to reconstruct the tree.)

ATGC5_parameters.txt - parameters fit the statistics of ATGC005

ATGC7_parameters.txt

ATGC8_parameters.txt

ATGC9_parameters.txt

ATGC32_parameters.txt

Average_25_parameters.txt - parameters fit the average of the participating ATGC, with 25 leaves

The next examples are parameter files that simulate one edge and test the ability of the tested methods to estimate the length:

one_edge_pure_cogs_jump_parameters.txt

one_edge_pure_cogs_parameters.txt

one_edge_pure_jump_parameters.txt

one_edge_pure_jump_parameters.txt

Output: printed on the screen.

2. Indel_real_data_creates_trees.py
   
This program creates a tree for a given ATGC number for each of the tested methods (DCJ, GC, GCU, or CGA)

For this, it opens one of the following files:

ATGC5reduced.csv

ATGC7reduced.csv

ATGC8reduced.csv

ATGC9reduced.csv

ATGC32reduced.csv

Example for the command line: 

py Indel_real_data_creates_trees.py  5

Create trees for ATGC005 using all the tested methods.

Another example:

py Indel_real_data_creates_trees.py  --No_DCJ 5

Similar to the previous command, but without DCJ in case that UNIMOG and Java are not available

Output: print the tree graph on the screen and create a Newick file for each method.

Newick files that were created by this function:

Files for the DCJ method:

DCJ_5_tree_data.txt

DCJ_7_tree_data.txt

DCJ_8_tree_data.txt

DCJ_9_tree_data.txt

DCJ_32_tree_data.txt

Files for the GC method:

GC_5_tree_data.txt

GC_7_tree_data.txt

GC_8_tree_data.txt

GC_9_tree_data.txt

GC_32_tree_data.txt

Files for the GCU method:

GCU_5_tree_data.txt

GCU_7_tree_data.txt

GCU_8_tree_data.txt

GCU_9_tree_data.txt

GCU_32_tree_data.txt

Files for the CGA method:

CGA_5_tree_data.txt

CGA_7_tree_data.txt

CGA_8_tree_data.txt

CGA_9_tree_data.txt

CGA_32_tree_data.txt

Newick files downloaded from the ATGCs site:

ATGC_5_tree_data.txt

ATGC_7_tree_data.txt

ATGC_8_tree_data.txt

ATGC_9_tree_data.txt

ATGC_32_tree_data.txt

3. Indel_Compare_Trees_Real_Data.py

This program measures, for a given ATGC number, the Normalized Robinson Foulds Distance (NRFD) between the trees created by the different methods and between the ATGC tree downloaded from the ATGCs site.

It uses the Newick files created by the previous program.

Example for the command line: 

py Indel_Compare_Trees_Real_Data.py  5

It compares the trees for ATGC005 for all the tested methods.

Another example:

py Indel_real_data_creates_trees.py  --No_DCJ 5

Similar to the previous command, but without DCJ in case that UNIMOG and Java are not available

Output: The results are printed on the screen.
5. Indel_real_data_trees_statistics.py

This program gives statistics for a given ATGC number, including statistics about the GCU and CGA trees.  It uses the Newick files created by the program in 2.

Example for the command line: 

py Indel_real_data_trees_statistics.py  5

It gives statistics about the ATGC005 family.

Output: The results are printed on the screen.

Tree graph files that can be found here:

dcj5.png

dcj7.png

dcj8.png

dcj9.png

dcj32.png

gc5.png

gc7.png

gc8.png

gc9.png

gc32.png

gcu5.png

gcu7.png

gcu8.png

gcu9.png

gcu32.png

cga5.png

cga7.png

cga8.png

cga9.png

cga32.png

atgc5.png

atgc7.png

atgc8.png

atgc9.png

atgc32.png















