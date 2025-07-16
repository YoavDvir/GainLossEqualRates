import math
import random
import statistics
import warnings
import pandas as pd

from numpy.ma.core import zeros
from scipy.optimize import fsolve

from itertools import product,filterfalse,starmap,islice,count
from timeit import default_timer as timer
import os

import numpy as np
from ete3 import Tree
from statistics import mean,stdev
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix

import argparse

def get_args():
    ind_sim = IndelSimulations()
    parser = argparse.ArgumentParser(description='perform simulation under gain-loss model')
    parser.add_argument('parameters_file',
                        metavar ='parameters_file',
                        type = str,
                        nargs = 1,
                        help ='path to parameters file')
    parser.add_argument(dest='sim',
                        action='store_const',
                        const= ind_sim.simulations,
                        help='run the program for the specified parameters file')
    args = parser.parse_args()
    return args.sim(args.parameters_file[0])
#Simulations of the gain-loss model under the equal rates assumption
#using Chinese restaurant process and jumps
#methods for measuring distance between genomes: cga,dcj

def only_first(l): #key function for sorting, according to the first element
    return l[0]

def only_second(l): #key function for sorting, according to the second element
    return l[1]

def cog_only(l): #taking only the cog names
    return list(map(only_first, l))

#
def set_of_pairs(gm): #create a set of pairs of consecutive genes
    l=len(gm)
    gm_set=set()
    for i in range(l-1):
        gm_set.add((gm[i], gm[i + 1]))
    return gm_set

def calc_dist_from_cga(t, cga):
    res= cga[0] - math.exp(-2 * t[0]) / (1 + t[0])
    return res

def create_match(pair0,pair1):
    #function for the lcs method: help to create a new list with an element for each pair of pairs that have the same first
    #element
    if pair0[0]==pair1[0]:
        return[pair0[0],pair0[1]+pair1[1],pair0[1],pair1[1],0,0]
    else:
        return []

def cr_index(lst): #zip to each element an index
    return list(zip(lst,islice(count(),len(lst))))

def print_dm(dm): #print distance matrix
    l=len(dm)
    for i in range(l):
        print(dm[i])

# distance functions
def dcj_dist(genome_pair): #dcj method (unimog package)
    err_flag = False
    leaf_genomes = list(genome_pair)
    try:
        os.remove('dcj_matrix.txt')
    except OSError:
        print("no file")

    try:
        os.remove('UNIMOG_file.txt')
    except OSError:
        print("no file")

    f_inp = open('UNIMOG_file.txt', 'w')
    leaves_number = 2
    for m in range(2):
        f_inp.write("a>" + str(m) + "\n")
        for COG in leaf_genomes[m]:
            f_inp.write(str(COG) + "\n")
    f_inp.close()
    os.system("java -jar UniMoG.jar -m=6 -d UNIMOG_file.txt >>dcj_matrix.txt")
    f = open("dcj_matrix.txt", "r")
    k = 0
    while k < 2:
        line = f.readline()
        lst = line.split()
        if len(lst) == 1:
            k += 1
    dm = []
    for m in range(leaves_number):
        dm.append([])
    for m in range(1, leaves_number):
        line = f.readline()
        line = line[10:]
        dm[m] = list(map(int, line.split()))
    for m in range(leaves_number):
        dm[m].append(0)
    return dm[1][0], err_flag

def gc_dist(genome_pair): #gc method
    err_flag = False
    g_set1 = genome_pair[0][1]
    lg_set1 = len(g_set1)
    g_set2 = genome_pair[1][1]
    lg_set2 = len(g_set2)
    intersect12 = g_set1.intersection(g_set2)
    l_intersect12 = len(intersect12)
    dist_sgc = -math.log(l_intersect12/min(lg_set1,lg_set2))
    return dist_sgc, err_flag

def gcu_dist(genome_pair): #gcu method
    err_flag = False
    g_set1 = genome_pair[0][1]
    lg_set1 = len(g_set1)
    g_set2 = genome_pair[1][1]
    lg_set2 = len(g_set2)
    intersect12 = g_set1.intersection(g_set2)
    l_intersect12 = len(intersect12)
    dist_cgc = -math.log(2*l_intersect12/(lg_set1 + lg_set2))
    return dist_cgc, err_flag

def cga_dist(genome_pair): #cga method

    err_flag = False
    g_set_pairs0 = genome_pair[0][2]
    lg_set_pairs0 = len(g_set_pairs0)

    g_set_pairs1 = genome_pair[1][2]
    lg_set_pairs1 = len(g_set_pairs1)

    intersect_pairs01 = g_set_pairs0.intersection(g_set_pairs1)
    l_intersect_pairs01 = len(intersect_pairs01)
    size_pairs_mean = (lg_set_pairs0 + lg_set_pairs1) / 2
    if l_intersect_pairs01 == 0:
        l_intersect_pairs01 = 1
    cga = l_intersect_pairs01 / size_pairs_mean

    with warnings.catch_warnings(record=True) as w:
        dist_cga = fsolve(calc_dist_from_cga, zeros(1), [cga])[0]
        if len(w) > 0:
            print("cga_error")
            dist_cga = 0
            err_flag = True
    return dist_cga, err_flag

def lcs_dist(genome_pair): #lcs method

    err_flag = False
    groups_index = list(map(cr_index, [genome_pair[0][0], genome_pair[1][0]]))
    prod_lst = list(product(groups_index[0], groups_index[1]))
    lst = list(filterfalse(lambda x: x == [], starmap(calc_dist_from_cga, prod_lst)))
    lst.sort(key=only_second)

    lst[0][4] = 1
    lst[0][5] = 0
    max_len_0up = 1
    max_len_1up = 0
    for i in range(1, len(lst)):
        for j in range(0, i):
            if lst[i][2] > lst[j][2] and lst[i][3] > lst[j][3]:
                if lst[i][4] < lst[j][4] + 1:
                    lst[i][4] = lst[j][4] + 1
                if lst[i][5] < lst[j][5]:
                    lst[i][5] = lst[j][5]
                if lst[i][2] == lst[j][2] + 1 and lst[i][3] == lst[j][3] + 1 and lst[i][5] < lst[j][5] + 1:
                    lst[i][5] = lst[j][5] + 1
        if lst[i][4] > max_len_0up:
            max_len_0up = lst[i][4]
        if lst[i][5] > max_len_1up:
            max_len_1up = lst[i][5]
    l1 = len(genome_pair[0])
    l2 = len(genome_pair[1])
    dist_lcs = -math.log(max_len_0up / (l1 * l2) ** 0.5)
    return dist_lcs, err_flag

class IndelSimulations:
    #program parameters
    #setting parameters
    trials = 100
    rt_size = 1800   #root size
    mean_indel = 0.03 #mean length of edge due to indels
    jmp = True # if True: with jumps
    mean_jmp = 0.075 # mean length of edge due to jumps
    chinese = True # if True: with nontrivial cogs
    theta = 717  # Chinese restaurant concentration parameter
    alpha = 0.965  # Chinese restaurant discount parameter
    w_dcj = True # True if with dcj method
    w_gc = True  # True if with gc method
    w_gcu = True  # True if with gcu method
    w_cga = True  # True if with cga method
    trees =  True # True if trees simulations. False if one-edge simulations
    leave_n = 25 #number of leaves of the tree
    std_indel = 0.03 #stdev of length of indel edges
    std_jmp = 0.155  # stdev of length of jmp edges
    # ATGC005
    """
    leaves_number = 21
    root_size = 1962
    theta = 581
    alph = 0.976
    mean_edge_lengths_indel = 0.023
    stdev_edge_lengths_indel = 0.056
    mean_edge_lengths_jump = 0.034
    stdev_edge_lengths_jump = 0.056
    """

    # ATGC007
    """
    leaves_number = 14
    root_size = 2333
    theta = 661
    alph = 0.949
    mean_edge_lengths_indel = 0.041
    stdev_edge_lengths_indel = 0.033
    mean_edge_lengths_jump = 0.027
    stdev_edge_lengths_jump = 0.048
    """

    # ATGC008
    """
    leaves_number = 15
    root_size = 1786
    theta = 746
    alph = 0.978
    mean_edge_lengths_indel = 0.029
    stdev_edge_lengths_indel = 0.039
    mean_edge_lengths_jump = 0.048
    stdev_edge_lengths_jump = 0.264
    """

    # ATGC009
    """
    leaves_number = 9
    root_size = 1916
    theta = 1089
    alph = 0.952
    mean_edge_lengths_indel = 0.051
    stdev_edge_lengths_indel = 0.033
    mean_edge_lengths_jump = 0.108
    stdev_edge_lengths_jump = 0.413
    """
    # ATGC032
    """
    leaves_number = 11
    root_size = 728
    theta = 704
    alph = 0.959
    mean_edge_lengths_indel = 0.010
    stdev_edge_lengths_indel = 0.014
    mean_edge_lengths_jump = 0.030
    stdev_edge_lengths_jump = 0.085
    """
    # Average
    """
    leaves_number = 25
    root_size = 1800
    theta = 717
    alph = 0.965
    mean_edge_lengths_indel = 0.03
    stdev_edge_lengths_indel = 0.03
    mean_edge_lengths_jump = 0.075
    stdev_edge_lengths_jump = 0.155
    """

    if trees:
        nrfd_max = 2 * leave_n - 6  # Maximal possible Robinson - Foulds distance
        shape_indel = (mean_indel / std_indel) ** 2
        scale_indel = mean_indel / shape_indel
        shape_jump = (mean_jmp / std_jmp) ** 2
        scale_jump = mean_jmp / shape_jump



    n0 = 0  # length of the genome at time t=0
    cogs_n = 0  # number of different cogs
    cogs = []  # Pool of all Cogs with repetitions
    MoreThanOne = []  # genes added to existing cogs
    OneMemberCogs = []  # Cogs with one gene
    leaves = [] # list of nodes which are leaves

    indel_jump_probabilities = [] # for each edge the proportion between indel edge length and total edge length
    indel_jump_probability = 0 # the proportion between indel edge length and total edge length for current edge
    rooted_tree_edge_lengths = [] #will contain the list of the random edge lengths for the rooted tree
    edge_lengths =[]    #will contain the list of the random edge lengths for the unrooted tree
    edge_length = []    #will contain a specific edge length

    random_tree = [] # structure: list of lists: each list contains:
                    # serial number,
                    # child 1 serial number(=0 if it is a leaf),
                    # child 2 serial number (=0 if it is a leaf),
                    # and parent serial number (=-1 if it is a root)
    genomes = [] #list of list: for each node (genome) the list of genes
    genomes_c = [] #list of list: for each node (genome) the list of cogs. the same shape as genomes
    node_serial_number = 0 # when the tree is created contains the new node serial number
    gene_serial = 0 # when the tree is created contains the new gene serial number
    leaf_genomes = [] # The output of the creation for the random tree
                        # and the input for the construction methods.
                        # list of list: for each leaf genome contains a list of genes
    leaf_genomes_c = [] # The output of the creation for the random tree
                        # and the input for the construction methods.
                        # list of list: for each leaf genome contains a list of cogs
    leave_sizes = [] #list of leaf sizes (number of genes)
    means_leave_size = [] #mean of leaf sizes (number of genes)
    stdev_leave_size = [] #stdev of leaf sizes (number of genes)
    mean_means_leave_size = 0 #mean of means leaf sizes (number of genes) for all repetitions
    stdev_means_leave_size = 0 #stdev of means leaf sizes (number of genes) for all repetitions
    mean_stdev_leave_size = 0 #mean of stdev leaf sizes (number of genes) for all repetitions
    stdev_stdev_leave_size = 0 #stdev of stdev leaf sizes (number of genes) for all repetitions
    methods = [] #method names
    mat_functions = [] #distance matrix for each tested method
    dist_functions = [] #distance between two genomes for each tested method
    l_functions = 0 #number of tested methods
    labels =[] #list of labels for the nodes

    constructed_trees = [] #list of constructed trees, one for each method (pure)
    constructed_trees_c = [] #list of constructed trees, one for each method (cogs)
    nrfd = [] #normalized RF distance for each method (pure)
    nrfd_c = [] #normalized RF distance for each method (cogs)
    nrfd12 = [] #normalized RF distance between methods (pure)
    nrfd12_c = [] #normalized RF distance between methods (cogs)
    times = []
    G0 = []
    G0_c = []
    G1 = []
    G1_c = []
    dist = []
    dist_c = []
    t = []
    t_c = []

    def create_random_tree(self):
        # creates the random tree
        self.random_tree = [[0, 0, 0, -1]]
        genome_0 = list(range(self.n0)) #root genome
        self.genomes = [genome_0]
        genome_c_0 = list(self.chinese_first_genome()) #chinese root genome
        self.genomes_c = [genome_c_0]
        self.gene_serial = self.n0
        self.leaves = [0]
        self.node_serial_number = 1
        current_number_of_leaves = 1
        number_of_edges = 2 * self.leave_n - 2
        # random length of edges
        indel_edges_lengths = np.random.gamma(self.shape_indel, scale = self.scale_indel, size=number_of_edges )
        if self.jmp:
            jump_edges_lengths = np.random.gamma(self.shape_jump, scale = self.scale_jump, size=number_of_edges )

        self.rooted_tree_edge_lengths = []
        self.edge_lengths = []
        if self.jmp:
            for i in range(number_of_edges):
                self.rooted_tree_edge_lengths.append(indel_edges_lengths[i] + jump_edges_lengths[i] / 2)
                self.indel_jump_probabilities.append(indel_edges_lengths[i] / self.rooted_tree_edge_lengths[i])
                self.edge_lengths.append(indel_edges_lengths[i] + jump_edges_lengths[i])

        else:
            for i in range(number_of_edges):
                self.rooted_tree_edge_lengths.append(indel_edges_lengths[i])
                self.edge_lengths.append(indel_edges_lengths[i])



        self.edge_lengths[0] += self.rooted_tree_edge_lengths[1]
        self.edge_lengths[1] += self.rooted_tree_edge_lengths[0]

        while current_number_of_leaves < self.leave_n:  # construct the random tree with leaves_number leaves
        # Eliminating parent from the leaf list
            loc = random.randrange(current_number_of_leaves)
            random_new_parent = self.leaves[loc]
            self.leaves.pop(loc)
            # creating two random children and adding them to the leaf list
            child_number = 0
            self.create_child(child_number,random_new_parent)
            child_number = 1
            self.create_child(child_number,random_new_parent)
            current_number_of_leaves = len(self.leaves)
        return

    def create_child(self,child_number,random_new_parent):
        parent_genome = list(self.genomes[random_new_parent])
        parent_genome_c = list(self.genomes_c[random_new_parent])
        self.edge_length = self.rooted_tree_edge_lengths[self.node_serial_number - 1]
        if self.jmp:
            self.indel_jump_probability = self.indel_jump_probabilities[self.node_serial_number - 1]
        next_genome , next_genome_c = self.create_next_genomes(parent_genome, parent_genome_c)
        self.genomes.append(list(next_genome))
        self.genomes_c.append(list(next_genome_c))
        self.random_tree[random_new_parent][1+child_number] = self.node_serial_number

        if self.node_serial_number < 3:
            self.random_tree.append([self.node_serial_number, 0, 0, 2 - child_number])
            self.leaves.append(self.node_serial_number)
            self.node_serial_number += 1
        else:
            self.random_tree.append([self.node_serial_number , 0, 0, random_new_parent])
            self.leaves.append(self.node_serial_number)
            self.node_serial_number +=1
        return

    #creating first Chinese genome
    def chinese_first_genome(self):
        self.cogs_n = 0               #number of different cogs
        genome_c = []               #Genome list of cogs with repetitions
        self.MoreThanOne= []         #genes added to existing cogs
        self.OneMemberCogs = []     #Cogs with one gene
        self.cogs = []
        for n in range(self.n0):
            p = random.random()
            n_cogs=len(self.cogs)        #Number of cogs with repetitions
            if p < (self.theta + self.cogs_n * self.alpha) / (self.theta + n_cogs):
                #first appearance of a cog
                self.cogs_n +=1
                self.cogs.append(self.cogs_n)
                self.OneMemberCogs.append(self.cogs_n)
                genome_c.append(self.cogs_n)
            else:
                #adding cog already exist
                if p > (n_cogs + self.theta - self.cogs_n * (1 - self.alpha)) / (self.theta + n_cogs):
                    cog = random.randrange(self.cogs_n)
                else:
                    cog = random.choice(self.MoreThanOne)
                genome_c.append(cog)
                self.MoreThanOne.append(cog)
                if cog in self.OneMemberCogs:
                    place=self.OneMemberCogs.index(cog)
                    self.OneMemberCogs.pop(place)
        return genome_c

    # create next genomes (regular and chinese)
    def create_next_genomes(self,genome, genome_c):
    #inputs:
    #genome - list of genes (pure)
    #genome_c - list of cogs
        next_genome = list(genome)
        next_genome_c = list(genome_c)
        lng = len(next_genome)
        t = np.random.exponential(scale = 1 / (2 * lng ), size = 1)[0]
        while t < self.edge_length:
            random_place = random.randrange(lng) #location of event
            p0 = random.random()
            if p0 < self.indel_jump_probability or not self.jmp : #indel (not a jump)
                p1 = random.random()
                if p1 < 0.5 : #equal rates of insertions and deletions
                    #insertion of gene
                    left_right = random.randint(0,1)
                    insertion_place = random_place + left_right
                    next_genome.insert(insertion_place, self.gene_serial)

                    #Chinese insertion
                    if self.chinese:
                        p2 = random.random()
                        if p2 < (self.theta + self.cogs_n * self.alpha)/(self.theta + self.gene_serial):
                            #new cog
                            self.cogs.append(self.cogs_n)
                            self.OneMemberCogs.append(self.cogs_n)
                            next_genome_c.insert(insertion_place, self.cogs_n)
                            self.cogs_n += 1
                        else:
                            #adding existing cog
                            if p2 > (self.gene_serial + self.theta - self.cogs_n * (1 - self.alpha)) / (self.theta + self.gene_serial):
                                #Choose one of the existing cogs. Each existing cog has an equal chance to be chosen
                                # in probability (1 - self.alph)) / (self.theta + self.gene_serial)
                                cog = random.randrange(self.cogs_n)
                            else:
                                # choose only among non-trivial cogs
                                cog = random.choice(self.MoreThanOne)

                            if cog in self.OneMemberCogs:
                                place=self.OneMemberCogs.index(cog)
                                self.OneMemberCogs.pop(place)

                            self.cogs.append(cog) #list of cogs with repetitions
                            self.MoreThanOne.append(cog)
                            next_genome_c.insert(insertion_place, cog)
                    self.gene_serial += 1
                else: #deletion
                    next_genome.pop(random_place)
                    if self.chinese:
                        next_genome_c.pop(random_place)
            else: #jump
                if self.jmp:
                    jumping_cog = next_genome[random_place]
                    if self.chinese:
                        jumping_cog_c = next_genome_c[random_place]
                    next_genome.pop(random_place)
                    if self.chinese:
                        next_genome_c.pop(random_place)
                    random_place2 = random.randrange(lng) #second location for the jump
                    next_genome.insert(random_place2, jumping_cog)
                    if self.chinese:
                        next_genome_c.insert(random_place2, jumping_cog_c)

            lng = len(next_genome)
            t += np.random.exponential(scale=1 / (2 * lng ), size=1)[0]

        return list(next_genome), list(next_genome_c)


#Constructing newick string and an ete unrooted tree given the rooted tree structure
    def newick(self):
        #inputs:
        # self.random_tree - the true tree structure
        # self.leaves - list of nodes numbers for leaf nodes
        # self.labels - list of nodes labels
        # self.leaves_number -number of leaves
        # outputs:
        # test_tree - true ete unrooted tree
        # potential output: 
        #true_dm - the true data matrix

        newick_tree = list(self.labels)
        root_subtree_serials = list(range(self.leave_n))
        new_node_serial = self.leave_n
        constructed_nodes = list(self.leaves)
        components_number = self.leave_n

        while components_number  > 3:
            pair_found = 0
            for m in range(components_number - 1):
                for n in range(m + 1 , components_number):
                    # check if the root_subtree_serials[m] and the root of the subtree 
                    # [root_subtree_serials[n]] have a common neighbor in the true tree
                    cna = constructed_nodes[root_subtree_serials[m]]
                    cnb = constructed_nodes[root_subtree_serials[n]]
                    parent = self.random_tree[cna][3]
                    child1 = self.random_tree[cna][1]
                    child2 = self.random_tree[cna][2]
                    neighbors1 = {parent}
                    if child1 > 0:
                        neighbors1.add(child1)
                    if child2 > 0:
                        neighbors1.add(child2)
                    parent = self.random_tree[cnb][3]
                    child1 = self.random_tree[cnb][1]
                    child2 = self.random_tree[cnb][2]

                    neighbors2 = {parent}
                    if child1 > 0:
                        neighbors2.add(child1)
                    if child2 > 0:
                        neighbors2.add(child2)
                    connect = neighbors1 & neighbors2
                    connect_list = list(connect)

                    if len(connect_list) == 1:
                        cnc = connect_list[0]
                        if not(cnc in constructed_nodes):
                            #two roots of subtrees have a common neighbor
                            pair_found = 1 
                            constructed_nodes.append(cnc)
                            newick_tree.append("(" + newick_tree[m] + "," + newick_tree[n]+")")
                            root_subtree_serials.append(new_node_serial)
                            
                            newick_tree.pop(n)
                            newick_tree.pop(m)
                            root_subtree_serials.pop(n)
                            root_subtree_serials.pop(m)
                            new_node_serial += 1
                            break
                if pair_found == 1:
                    break
            if pair_found == 0:
                print("error")
                exit()
            components_number = len(root_subtree_serials)
        newick_tree_string = "(" + newick_tree[0] + "," + newick_tree[1] + "," + newick_tree[2] + ");"

        test_tree =Tree(newick_tree_string)
        return test_tree

    #construct a tree from a distance matrix
    def dm_to_nj_tree(self , dm):
        m = DistanceMatrix(list(self.labels), list(dm))
        constructor = DistanceTreeConstructor()
        nj_tree = constructor.nj(m)
        tree_data = "treedata.txt"
        Phylo.write(nj_tree, tree_data, "newick")
        f = open(tree_data, "r")
        newick_tree = f.read()
        constructed_tree = Tree(newick_tree, format=1)
        return constructed_tree

    #methods to get distance matrices from a list of genome leaves
    # methods: djc, gcu, gc, cga,

    #compute djc distance matrix from a list of leaf genomes
    def djc_mat(self,local_leaf_genomes):
        try:
            os.remove('dcj_matrix.txt')
        except OSError:
            print("no file")
        try:
            os.remove('UNIMOG_file.txt')
        except OSError:
            print("no file")
        f_inp = open('UNIMOG_file.txt', 'w')
        for m in range(self.leave_n):
            f_inp.write("a>" + str(m) + "\n")
            for COG in local_leaf_genomes[m][0]:
                f_inp.write(str(COG) + "\n")
        f_inp.close()
        os.system('java -jar UniMoG.jar -m=6 -d UNIMOG_file.txt >>dcj_matrix.txt')
        f = open('dcj_matrix.txt', 'r')
        k = 0
        while k < 3:
            line = f.readline()
            lst = line.split()
            if len(lst) == 1:
                k += 1
        dm= []
        for m in range(self.leave_n):
            dm.append([])
        for m in range(1, self.leave_n):
            line = f.readline()
            line = line[10:]
            dm[m] = list(map(int, line.split()))
        for m in range(self.leave_n):
            dm[m].append(0)
        return dm
    # compute gc distance matrix from a list of leaf genomes
    def gc_mat(self,g):
        dm = []
        for m in range(self.leave_n):
            dm.append([])

        for m in range(self.leave_n - 1):
            for n in range(m + 1, self.leave_n):
                pairs_content, err_flag = gc_dist([g[m], g[n]])
                if err_flag:
                    print('pairs_content error')
                    return None
                dm[n].append(pairs_content)
        for m in range(self.leave_n):
            dm[m].append(0)
        return dm

    # compute gcu distance matrix from a list of leaf genomes
    def gcu_mat(self,g):
        dm = []
        for m in range(self.leave_n):
            dm.append([])

        for m in range(self.leave_n - 1):
            for n in range(m + 1, self.leave_n):
                pairs_content, err_flag = gcu_dist([g[m], g[n]])
                if err_flag:
                    print('pairs_content error')
                    return None
                dm[n].append(pairs_content)
        for m in range(self.leave_n):
            dm[m].append(0)
        return dm

    # compute cga distance matrix from a list of leaf genomes
    def cga_mat(self,g):
        dm = []
        for m in range(self.leave_n):
            dm.append([])

        for m in range(self.leave_n - 1):
            for n in range(m + 1, self.leave_n):
                pairs_content, err_flag = cga_dist([g[m], g[n]])
                if err_flag:
                    print('pairs_content error')
                    return None
                dm[n].append(pairs_content)
        for m in range(self.leave_n):
            dm[m].append(0)
        return dm

    def distance_for_all_method(self):
        #simulate one edge and measure distance by all method
        #creates the root genome (pure)
        err_flag = False
        genome_0 = list(range(self.n0))  # root genome
        genome_s_0 = set(genome_0)
        genome_p_0 = set_of_pairs(genome_0)
        self.G0 = [genome_0, genome_s_0, genome_p_0] #keep the sequence of genes, set of genes and set of consecutive pairs

        # creates the root genome (cogs)
        if self.chinese:
            genome_c_0 = list(self.chinese_first_genome())  # chinese root genome
            genome_s_c_0 =set(genome_c_0)
            genome_p_c_0 = set_of_pairs(genome_c_0)
        else:
            genome_c_0 = []
            genome_s_c_0 = set()
            genome_p_c_0 = set()
        self.G0_c = [genome_c_0, genome_s_c_0, genome_p_c_0]

        #create the genomes after the simulated edge
        self.gene_serial = self.n0
        next_genome , next_genome_c = self.create_next_genomes(genome_0, genome_c_0)

        #pure
        genome_1 = list(next_genome)
        genome_s_1 = set(genome_1)
        genome_p_1 = set_of_pairs(genome_1)
        self.G1 = [genome_1, genome_s_1, genome_p_1]

        #cogs
        if self.chinese:
            genome_c_1 = list(next_genome_c)
            genome_s_c_1 = set(genome_c_1)
            genome_p_c_1 = set_of_pairs(genome_c_1)
            self.G1_c = [genome_c_1, genome_s_c_1, genome_p_c_1]

        for f in range(self.l_functions): #measure distances (pure and cogs) for each tested method
            start = timer()
            d, err_flag= self.dist_functions[f]([self.G0,self.G1])
            end = timer()
            dt = end - start
            self.t[f].append(dt)
            self.dist[f].append(d)
            if err_flag:
                print("error")
                return err_flag
            if self.chinese:
                start = timer()
                d_c, err_flag = self.dist_functions[f]([self.G0_c,self.G1_c])
                end = timer()
                dt_c = end - start
                self.t_c[f].append(dt_c)
                self.dist_c[f].append(d_c)
                if err_flag:
                    print("error")
                    return err_flag
        return err_flag

    def rfd_for_all_method(self):
        #input:
        #self.leaf_genomes - List of lists. For each leaf genome - list of genes (pure)
        #self.leaf_genomes_c - List of lists. For each leaf genome - list of genes (cogs)
        #self.random_tree (as input for self.newick)
        #output:
        #nrfd[f]
        #creates distance matrices by all method than compute trees that are compared to the true tree by nrfd
        #nrfd12[f][f2] - nfr distance between methods
        test_tree = self.newick() #get ete tree from self.random_tree
        err_flag = False
        genomes = list(self.leaf_genomes)
        if self.chinese:
            genomes_c = list(self.leaf_genomes_c)
        l_genomes =len(genomes)
        g = [] #will contain the leaf genomes data (pure)
        g_c = [] #will contain the leaf genomes data (cogs)
        for i in range(l_genomes): #for each leaf
            genome = list(genomes[i])
            genome_set = set(genomes[i])
            genome_pair_set = set_of_pairs(genome)
            g.append([genome, genome_set, genome_pair_set])
            if self.chinese:
                genome_c = list(genomes_c[i])
                genome_set_c = set(genomes_c[i])
                genome_pair_set_c = set_of_pairs(genome_c)
                g_c.append([genome_c, genome_set_c, genome_pair_set_c])
        self.constructed_trees = []
        self.constructed_trees_c = []
        for f in range(self.l_functions): #creates a distance matrix for each method (pure and cogs)
            start = timer()
            dm = self.mat_functions[f](g) #creates distance matrix (pure}
            end = timer()
            dt = end - start
            if dm is None:
                err_flag = True
                return err_flag
            self.t[f].append(dt)
            if self.chinese:
                start = timer()
                dm_c = self.mat_functions[f](g_c) # creates distance matrix (cogs)
                end = timer()
                dt = end - start
                if dm_c is None:
                    err_flag = True
                    return err_flag
                self.t_c[f].append(dt)
            #create ete trees from distance matrices
            result_nj_tree = self.dm_to_nj_tree(dm)
            self.constructed_trees.append(result_nj_tree)
            d = test_tree.robinson_foulds(result_nj_tree, unrooted_trees=True)[0]
            self.nrfd[f].append(d / self.nrfd_max)
            if self.chinese:
                result_nj_tree_c = self.dm_to_nj_tree(dm_c)
                self.constructed_trees_c.append(result_nj_tree_c)
                d_c = test_tree.robinson_foulds(result_nj_tree_c, unrooted_trees=True)[0]
                self.nrfd_c[f].append(d_c / self.nrfd_max)
            for f2 in range(f): #compare trees of different methods
                d = result_nj_tree.robinson_foulds(self.constructed_trees[f2], unrooted_trees=True)[0]
                self.nrfd12[f][f2].append(d / self.nrfd_max)
                if self.chinese:
                    d_c = result_nj_tree_c.robinson_foulds(self.constructed_trees_c[f2], unrooted_trees=True)[0]
                    self.nrfd12_c[f][f2].append(d_c / self.nrfd_max)
        return err_flag

    #perform repetitions of simulated edge for a given set of parameters and set of methods
    def repetitions_edge(self):

        self.dist = []
        self.dist_c = []

        self.t = []
        self.t_c = []



        dist_m = []
        dist_s = []

        dist_c_m = []
        dist_c_s = []

        t_m =[]
        t_m_c = []

        for f in range(self.l_functions):
            self.dist.append([])
            if self.chinese:
                self.dist_c.append([])
            self.t.append([])
            if self.chinese:
                self.t_c.append([])

            self.times.append([])

        print("self.trials",self.trials)
        error_number = 0
        for i in range(self.trials):
            print("trial",i)


            if self.distance_for_all_method():
                print("error")
                error_number += 1

        print("error_number", error_number)
        for f in range(self.l_functions):
            dist_mean = round(statistics.mean(self.dist[f]),6)
            dist_stdev = round(statistics.stdev(self.dist[f])/dist_mean,6)
            dist_m.append(str(dist_mean))
            dist_s.append(str(dist_stdev))

            if self.chinese:
                dist_c_mean = round(statistics.mean(self.dist_c[f]),6)
                dist_c_stdev = round(statistics.stdev(self.dist_c[f])/dist_c_mean,6)
                dist_c_m.append(str(dist_c_mean))
                dist_c_s.append(str(dist_c_stdev))

            t_mean = statistics.mean(self.t[f])
            t_m.append(t_mean)
            if self.chinese:
                t_mean = statistics.mean(self.t_c[f])
                t_m_c.append(t_mean)
        labels = ["mean distance (pure)", "stdev distance (pure)", "mean distance (cogs)", "stdev distance (cogs)",
                  "mean time (pure)","mean time (cogs)"]
        results = [dist_m, dist_s, dist_c_m, dist_c_s, t_m, t_m_c]
        return [results, labels]

    #perform repetitions of simulated random trees for a given set of parameters and set of methods
    def repetitions_trees(self):

        self.nrfd = []
        self.nrfd_c = []

        self.nrfd12 = []
        self.nrfd12_c = []

        t_m = []
        t_c_m =[]
        self.times = []

        nrfd_m = []
        nrfd_s = []
        nrfd12_m = []

        nrfd_c_m = []
        nrfd_c_s = []
        nrfd12_c_m = []

        nrfd12_tre_m = []
        nrfd12_tre_c_m = []

        for f in range(self.l_functions):
            self.nrfd.append([])
            self.nrfd12.append([])

            self.nrfd_c.append([])
            self.nrfd12_c.append([])

            self.t.append([])
            self.t_c.append([])

            nrfd12_m.append([])
            nrfd12_c_m.append([])

            nrfd12_tre_m.append([])
            nrfd12_tre_c_m.append([])

            for f2 in range(f):
                self.nrfd12[f].append([])
                self.nrfd12_c[f].append([])

        rf_matrix = []
        mean_rf_matrix = []
        sd_rf_matrix = []
        for f1 in range(4 * self.l_functions - 1):
            rf_matrix.append([])
            mean_rf_matrix.append([])
            sd_rf_matrix.append([])

        self.labels = []
        for m in range(self.leave_n):
            self.labels.append("L" + str(m))

        for f1 in range(4 * self.l_functions - 1):
            for f2 in range(f1 , 4 * self.l_functions - 1):
                rf_matrix[f2].append([])
                mean_rf_matrix[f2].append([])
                sd_rf_matrix[f2].append([])
        
        self.nrfd_max = 2 * self.leave_n - 6  # Maximal possible Robinson - Foulds distance
        self.shape_indel = (self.mean_indel / self.std_indel) ** 2
        self.scale_indel = self.mean_indel / self.shape_indel
        self.shape_jump = (self.mean_jmp / self.std_jmp) ** 2
        self.scale_jump = self.mean_jmp / self.shape_jump
        
        print("leaves number", self.leave_n)
        error_number = 0
        #print("self.trials", self.trials)
        for i in range(self.trials):
            #print("trial",i)
            self.constructed_trees = []
            self.constructed_trees_c = []

            self.create_random_tree()
            self.leaf_genomes = []
            self.leaf_genomes_c = []
            self.leave_sizes =  []

            for m in self.leaves:
                self.leaf_genomes.append(self.genomes[m])
                if self.chinese:
                    self.leaf_genomes_c.append(self.genomes_c[m])
                self.leave_sizes.append(len(self.genomes[m]))
            self.means_leave_size.append(np.mean(self.leave_sizes))
            self.stdev_leave_size.append(np.std(self.leave_sizes))


            if self.rfd_for_all_method():
                print("error")
                error_number += 1

        print("error_number", error_number)

        for f in range(self.l_functions):
            nrfd_m.append(statistics.mean(self.nrfd[f]))
            nrfd_s.append(statistics.stdev(self.nrfd[f]))
            nrfd_c_m.append(statistics.mean(self.nrfd_c[f]))
            nrfd_c_s.append(statistics.stdev(self.nrfd_c[f]))
            t_m.append(statistics.mean(self.t[f]))
            t_c_m.append(statistics.mean(self.t_c[f]))
            for f2 in range(f):
                nrfd12_m[f].append(statistics.mean(self.nrfd12[f][f2]))
                nrfd12_c_m[f].append(statistics.mean(self.nrfd12_c[f][f2]))

        self.mean_means_leave_size = mean(self.means_leave_size)
        self.stdev_means_leave_size = stdev(self.means_leave_size)
        self.mean_stdev_leave_size = mean(self.stdev_leave_size)
        self.stdev_stdev_leave_size = stdev(self.stdev_leave_size)

        labels = ["mean nrfd (pure)", "stdev nrfd (pure)", "mean nrfd (cogs)",
                  "stdev nrfd (cogs)","mean time (pure)", "mean time (cogs)"]

        results = [nrfd_m, nrfd_s, nrfd_c_m, nrfd_c_s, t_m, t_c_m]
        df=pd.DataFrame(nrfd12_m[1:self.l_functions],columns=self.methods[0:self.l_functions-1],index=self.methods[1:self.l_functions])
        print("comparing nrfd between trees of different methods:")
        print(df)
        return [results, labels]

    def simulations(self,file_name): #main program for simulations of "one edge" or "trees"
        total_time_start = timer()
        # program parameters
        #df = pd.read_csv('one_edge_pure_cogs_jump_parameters.txt')
        #df = pd.read_csv('one_edge_pure_jump_parameters.txt')
        #df = pd.read_csv('one_edge_pure_cogs_parameters.txt')
        #df = pd.read_csv('ATGC9_parameters.txt')

        df = pd.read_csv(file_name)
        print("parameters set:")
        print(df.to_string())
        # setting parameters
        self.trials = df.loc[0, "trials"]  # number of repetitions
        self.rt_size = df.loc[0, "rt_size"]  # root size
        self.mean_indel = df.loc[0, "mean_indel"]  # mean length of edge due to indels
        self.jmp = df.loc[0, "jmp"]  # if True: with jumps
        self.mean_jmp = df.loc[0, "mean_jmp"]  # mean length of edge due to jumps
        self.chinese = df.loc[0, "chinese"]  # if True: with nontrivial cogs
        self.theta = df.loc[0, "theta"]  # Chinese restaurant concentration parameter
        self.alpha = df.loc[0, "alpha"]  # Chinese restaurant discount parameter
        self.l_functions = 0
        self.w_dcj = df.loc[0, "w_dcj"]  # True if with dcj method
        if self.w_dcj:
            self.methods.append("DCJ")
            self.dist_functions.append(dcj_dist)
            self.mat_functions.append(self.djc_mat)
            self.l_functions =+1

        self.w_gc = df.loc[0, "w_gc"]  # True if with gc method
        if self.w_gc:
            self.methods.append("GC")
            self.dist_functions.append(gc_dist)
            self.mat_functions.append(self.gc_mat)
            self.l_functions += 1
        self.w_gcu = df.loc[0, "w_gcu"]  # True if with gcu method
        if self.w_gcu:
            self.methods.append("GCU")
            self.dist_functions.append(gcu_dist)
            self.mat_functions.append(self.gcu_mat)
            self.l_functions += 1
        self.w_cga = df.loc[0, "w_cga"]  # True if with cga method
        if self.w_cga:
            self.methods.append("CGA")
            self.dist_functions.append(cga_dist)
            self.mat_functions.append(self.cga_mat)
            self.l_functions += 1
        self.trees = df.loc[0, "trees"]  # True if trees simulations. False if one-edge simulations
        self.leave_n = df.loc[0, "leave_n"]  # number of leaves of the tree
        self.std_indel = df.loc[0, "std_indel"]  # stdev of length of indel edges
        self.std_jmp = df.loc[0, "std_jmp"]  # stdev of length of jmp edges

        self.n0 = self.rt_size

        #for i in range(len(self.mean_edge_lengths)):
        #    mel = self.mean_edge_length
        #    indel_jump_prob = self.mean_edge_prob[i]

        if self.trees:
            results, labels = self.repetitions_trees()
        else:
            if self.jmp:

                self.edge_length = self.mean_indel + self.mean_jmp / 2
                self.indel_jump_probability = self.mean_indel / self.edge_length
            else:
                self.edge_length = self.mean_indel

            results, labels = self.repetitions_edge()
        df=pd.DataFrame(results,columns=self.methods,index=labels)
        #df = pd.DataFrame(results)
        print(df)
        total_time_end = timer()
        print("Total time taken", total_time_end - total_time_start)
        return results, labels
get_args()
#ind_sim = IndelSimulations()
#ind_sim.dist_functions = [dcj_dist, gc_dist, gcu_dist, cga_dist]
#ind_sim.mat_functions = [ind_sim.djc_mat, ind_sim.gc_mat, ind_sim.gcu_mat, ind_sim.cga_mat]
#ind_sim.l_functions = len(ind_sim.mat_functions)
#ind_sim.mat_functions = [ind_sim.sgc_mat, ind_sim.cgc_mat]
#ind_sim.simulations()
