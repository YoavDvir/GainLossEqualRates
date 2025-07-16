
import pandas as pd
from ete3 import Tree
import argparse
#from Bio.Phylo import draw_ascii
def get_args():
    ind_real = IndelRealData()
    parser = argparse.ArgumentParser(description='compare trees from different methods.')
    parser.add_argument('ATGC_number',
                        metavar ='ATGC_number',
                        type = int,
                        nargs = 1,
                        help ='ATGC number')
    parser.add_argument('--No_DCJ', action='store_true', help='without DCJ')
    parser.add_argument(dest='compare',
                        action='store_const',
                        const= ind_real.real_data,
                        help='run the program for the specified ATGC number')
    args = parser.parse_args()
    return args.compare(args.ATGC_number[0],args.No_DCJ)
def print_dm(dm): #print distance matrix
    l=len(dm)
    for i in range(l):
        print(dm[i])
#read a tree for a method and atgc number
def read_tree(function_name, atgc_num): #read the tree that was created by the given method
    tree_data = function_name +"_"+ str(atgc_num) + '_tree_data.txt'
    f = open(tree_data, "r")
    newick_tr = f.read()
    constructed_tree = Tree(newick_tr, format=1)
    leaves_number = len(constructed_tree.get_leaves())
    return constructed_tree, leaves_number

class IndelRealData:

    leaves_number = 0
    nrfd_max = 0 # 2 * leaves_number - 6 # Maximal possible RF distance

    taxa_names = []
    node_serial_number = 0
    leaf_genomes = []

    functions = []
    functions_names = []
    l_functions = 0
    labels =[]

    edge_lengths =[]
    rooted_edge_lengths = []
    atgc_list = [(5, 22), (7, 14), (8, 15), (9, 9), (32, 11)]

    def real_data(self,atgc,no_dcj):
        print("No DCJ:", no_dcj)
        print("ATGC number: ", atgc)
        times = []
        if no_dcj:
            self.functions_names = ["ATGC","GC", "GCU", "CGA"]
            self.l_functions = 4
        else:
            self.functions_names = ["ATGC","DCJ", "GC", "GCU", "CGA"]
            self.l_functions = 5

        for f in range(self.l_functions):
            times.append([])
        nrfd_matrix = []

        for f1 in range(self.l_functions - 1):
            nrfd_matrix.append([])

        constructed_trees = []

        for f in range(self.l_functions):

            method_tree, leaves_n = read_tree(self.functions_names[f], atgc)
            constructed_trees.append(method_tree)
            
        self.labels = []
        for m in range(leaves_n):
            self.labels.append("L" + str(m))
        self.nrfd_max = 2 * leaves_n - 6

        for f1 in range(self.l_functions - 1):
            for f2 in range(f1 , self.l_functions - 1):
                d = constructed_trees[f1].robinson_foulds(constructed_trees[f2 + 1], unrooted_trees=True)[0]
                nrfd_matrix[f2].append(d / self.nrfd_max)


        print("nrfd matrix:")
        df=pd.DataFrame(nrfd_matrix,columns= self.functions_names[0:self.l_functions-1],index=self.functions_names[1:self.l_functions])
        print(df)


        return

get_args()

#ind_real = IndelRealData()

#ind_real.functions_names =["ATGC","DCJ","GC","GCU","CGA"]

#ind_real.l_functions = len(ind_real.functions_names)

#ind_real.real_data(self.atgc)
