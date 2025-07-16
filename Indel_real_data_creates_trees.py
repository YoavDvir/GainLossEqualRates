import math
from scipy.optimize import fsolve
from itertools import groupby
import os
import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio.Phylo import draw_ascii

#
import argparse

def get_args():
    ind_real = IndelRealData()
    parser = argparse.ArgumentParser(description='reconstruct trees for all methods given ATGC number')
    parser.add_argument('ATGC_number',
                        metavar ='ATGC_number',
                        type = int,
                        nargs = 1,
                        help ='ATGC number')
    parser.add_argument('--No_DCJ', action='store_true', help='without DCJ')
    parser.add_argument(dest='recon',
                        action='store_const',
                        const= ind_real.real_data,
                        help='run the program for the specified ATGC number')
    args = parser.parse_args()
    return args.recon(args.ATGC_number[0],args.No_DCJ)

def set_of_pairs(gm):# for CGA method, creates a set of consecutive pairs from gm
    l=len(gm)
    gm_set=set()
    for i in range(l-1):
        gm_set.add((gm[i], gm[i + 1]))
    return gm_set

def calc_dist_from_cga(t, cga): #given cga, the common number of consecutive pairs between two genomes
    #it creates an equation for solving the distance between two genomes under the model assumptions
    res= cga - math.exp(-2 * t[0]) / (1 + t[0])
    return res

def only_second(lst):
    return lst[1]

def only_first(lst):
    return lst[0]

def cog_only(lst): #take from the atgc file only the cogs column
    return list(map(only_first, lst))

def cga_distance(gm1,gm2): # measure the cga distance between two genomes

    gm_set_pairs1 = set_of_pairs(gm1)
    l_gm_set_pairs1 = len(gm_set_pairs1)

    gm_set_pairs2 = set_of_pairs(gm2)
    l_gm_set_pairs2 = len(gm_set_pairs2)

    intersect_pairs12 = gm_set_pairs1.intersection(gm_set_pairs2)
    l_intersect_pairs12 = len(intersect_pairs12)
    size_pairs_mean = (l_gm_set_pairs1 + l_gm_set_pairs2) / 2 - 1
    if l_intersect_pairs12 == 0:
        l_intersect_pairs12 = 1
    cga = l_intersect_pairs12 / size_pairs_mean
    dist_cga = fsolve(calc_dist_from_cga, 1, cga)[0]
    return dist_cga

class IndelRealData:

    leaves_number = 0
    nrfd_max = 2 * leaves_number - 6  # Maximal possible RF distance
    with_unimog = var = True

    taxa_names = []
    node_serial_number = 0
    leaf_genomes = []

    labels =[]

    edge_lengths =[]
    rooted_edge_lengths = []

    dm = []
    local_leaf_genomes = []
    local_leaf_genomes_sq = []
    cgc_stat = []
    cga_stat = []
    functions =[]
    functions_names =[]
    l_functions = len(functions_names)


    def read_atgc(self, atgc_num):
        #Reads the proper atgc file and fill the data to taxa_names, cogs_lists and n_leaves
        atgc_file = "ATGC" + str(atgc_num) + "reduced.csv"
        df = pd.read_csv(atgc_file, header=None)
        atgc_f = df.to_numpy()
        groups = []
        self.taxa_names = []

        atgc_f = sorted(atgc_f, key=only_second)

        for ky, gm in groupby(atgc_f, only_second):
            groups.append(list(gm))  # Store group iterator as a list
            self.taxa_names.append(ky)
        self.leaf_genomes = list(map(cog_only, groups))
        self.leaves_number = len(self.taxa_names)
        self.nrfd_max =  2 * self.leaves_number - 6
        self.labels = []
        for i in range(self.leaves_number):
            print(self.taxa_names[i])
            self.labels.append("L" + str(i))
        return

    #construct a nj tree (neighbor joining) from a distance matrix
    def dm_to_tree(self , dm , function_name , atgc_num):
        m = DistanceMatrix(list(self.taxa_names), list(dm))
        constructor = DistanceTreeConstructor()
        nj_tree = constructor.nj(m)
        print()
        print("reconstructed tree for method " + function_name + ":")
        print()
        draw_ascii(nj_tree)
        print()
        tree_data = function_name +"_"+ str(atgc_num) + "_tree_data.txt"
        Phylo.write(nj_tree, tree_data, "newick")
        return
    def dcj_mat(self):
        try:
            os.remove('dcj_matrix.txt')
        except OSError:
            print("no file")
        try:
            os.remove('UNIMOG_file.txt')
        except OSError:
            print("no file")
        f_inp = open('UNIMOG_file.txt', 'w')
        for m in range(self.leaves_number):
            f_inp.write("a>" + str(m) + "\n")
            for COG in self.local_leaf_genomes[m]:
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
        for m in range(self.leaves_number):
            dm.append([])
        for m in range(1, self.leaves_number):
            line = f.readline()
            line = line[10:]
            dm[m] = list(map(int, line.split()))
        for m in range(self.leaves_number):
            dm[m].append(0)
        return dm
    
    #compute distance matrix from a list of leaf genomes using the gene content method
    def gc_mat(self):
        dm = []
        cogsets = []
        for m in range(self.leaves_number):
            dm.append([])
            cogsets.append(set(list(self.local_leaf_genomes[m])))

        for m in range(self.leaves_number - 1):
            le0 = len(cogsets[m])
            for n in range(m + 1, self.leaves_number):
                le1 = len(cogsets[n])
                mn_intersection = cogsets[m].intersection(cogsets[n])
                mn_intersection_size = len(mn_intersection)
                if mn_intersection_size == 0:
                    mn_intersection_size = 1
                gene_content = -np.log(mn_intersection_size / min(le0, le1))
                dm[n].append(gene_content)
        for m in range(self.leaves_number):
            dm[m].append(0)

        return dm

    # compute distance matrix from a list of leaf genomes using the consistent probability gene content method
    def gcu_mat(self):
        dm = []
        cogsets = []
        for m in range(self.leaves_number):
            dm.append([])
            cogsets.append(set(list(self.local_leaf_genomes[m])))

        for m in range(self.leaves_number - 1):
            le0 = len(cogsets[m])
            for n in range(m + 1, self.leaves_number):
                le1 = len(cogsets[n])
                mn_intersection = cogsets[m].intersection(cogsets[n])
                mn_intersection_size = len(mn_intersection)
                if mn_intersection_size == 0:
                    mn_intersection_size = 1
                gene_content = -np.log(2 * mn_intersection_size / (le0 + le1))
                self.cgc_stat.append(gene_content)
                dm[n].append(gene_content)
        for m in range(self.leaves_number):
            dm[m].append(0)
        return dm


    # compute cga distance matrix from a list of leaf genomes
    def cga_mat(self):
        dm = []
        for m in range(self.leaves_number):
            dm.append([])

        for m in range(self.leaves_number - 1):
            for n in range(m + 1, self.leaves_number):
                pairs_content = cga_distance(self.local_leaf_genomes[m], self.local_leaf_genomes[n])
                self.cga_stat.append(pairs_content)
                dm[n].append(pairs_content)
        for m in range(self.leaves_number):
            dm[m].append(0)
        return dm

    #reconstruct a tree for a set of methods
    def real_data(self,atgc_num,no_dcj):
        print("ATGC", atgc_num)
        if no_dcj:
            self.functions = [self.gc_mat, self.gcu_mat, self.cga_mat]
            self.functions_names = ["GC", "GCU", "CGA"]
            self.l_functions = 3
        else:
            self.functions = [self.dcj_mat, self.gc_mat, self.gcu_mat, self.cga_mat]
            self.functions_names = ["DCJ", "GC", "GCU", "CGA"]
            self.l_functions = 4

        times = []

        for f in range(self.l_functions):
            times.append([])
        rf_matrix = []

        for f1 in range(self.l_functions - 1):
            rf_matrix.append([])

        self.labels = []
        for m in range(self.leaves_number):
            self.labels.append("L" + str(m))

        self.leaf_genomes = []

        self.read_atgc(atgc_num)

        self.local_leaf_genomes = list(self.leaf_genomes)

        mats = []
        for f in range(self.l_functions):
            self.dm = self.functions[f]()
            self.dm_to_tree(self.dm, self.functions_names[f], atgc_num)
            mats.append(list(self.dm))
        return



get_args()
#ind_real = IndelRealData()
#ind_real.functions = [ind_real.dcj_mat,ind_real.gc_mat, ind_real.gcu_mat, ind_real.cga_mat]
#ind_real.functions_names =["DCJ","GC","GCU","CGA"]
#ind_real.l_functions = len(ind_real.functions)
#ind_real.real_data()
