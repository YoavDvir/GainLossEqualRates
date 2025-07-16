from statistics import mean, stdev
import numpy as np
import argparse

def get_args():
    ind_stat = IndelRealDataStatistics()
    parser = argparse.ArgumentParser(description='collect statistics about gcu and cga trees.')
    parser.add_argument('ATGC_number',
                        metavar ='ATGC_number',
                        type = list,
                        nargs = 1,
                        help ='ATGC number')
    parser.add_argument(dest='stat',
                        action='store_const',
                        const= ind_stat.statistics,
                        help='run the program for the specified ATGC number')
    args = parser.parse_args()
    return args.stat(args.ATGC_number[0])
#read a tree for a method and atgc number
def read_tree_newick(function_name, atgc_num):
    treedata =  function_name + "_" + str(atgc_num) + "_tree_data.txt"
    f = open(treedata, "r")
    newick = f.read()
    return newick

def edge_length_list(newick):
    #fetch the edge lengths of the tree from the tree newick
    temp = newick
    length_list = []
    while len(temp) > 0:
        index = temp.find(":")
        if index == -1:
            break
        if temp[index + 1] == "-":
            numb = temp[index + 1: index + 9]
            temp = temp[ index + 8:]
        else:
            numb = temp[index + 1: index + 8]
            temp = temp[ index + 7:]
        fl_numb = float( numb )
        length_list.append(float(fl_numb))
    return length_list

class IndelRealDataStatistics:

    leaves_number = 0 # number of leaves in the tree
    nrfd_max = 0 # 2 * leaves_number - 6  # Maximal possible Robinson - Foulds distance
    with_unimog = var = True # using unimog or not

    taxa_names = [] # taxa names
    node_serial_number = 0 # next new node serial number
    leaf_genomes = [] # genomes of the leaves

    functions_names = [] # methods for the tree construction
    l_functions = 0 # number of methods for the tree construction

    ATGCs = [(5, 22), (7, 14), (8, 15), (9, 9), (32, 11)] # for each atgc: serial number, number of genomes
    tot_lst =[[],[]]

    def real_data(self,atgc):
        #print mean and stdev of edge lengths for each tree for a given atgc
        print("ATGC number:",atgc[0])
        for f in range(self.l_functions): #run for each method
            newick = read_tree_newick(self.functions_names[f],atgc[0]) #fetch the tree newick given the method and atgc number
            lst = edge_length_list(newick) #fetch the edge lengths of the tree
            self.tot_lst[f] += lst #create a list of lists of edge lengths for each tree
            #print("f",f,"atgc",atgc,"mean length", mean(lst), "stdev of length",stdev(lst)) # print statistics for each tree
        return
    def print_res(self):
        #print mean and stdev for cgu,cga and correlation between cgu and cga-cgu trees edge lengths
        # for each tree given atgc
        gcu_mean = mean(self.tot_lst[0])
        cga_mean = mean(self.tot_lst[1])
        gcu_stdev = stdev(self.tot_lst[0])
        cga_stdev = stdev(self.tot_lst[1])
        gcu_array = np.array(self.tot_lst[0])
        cga_array = np.array(self.tot_lst[1])
        diff_array = cga_array - gcu_array
        corr = np.corrcoef(gcu_array, diff_array)

        print("mean GCU length:", gcu_mean, "stdev:", gcu_stdev)
        print("mean CGA length:",cga_mean,"stdev:", cga_stdev)
        print("mean difference:", cga_mean - gcu_mean, "stdev:", (cga_stdev ** 2 - gcu_stdev ** 2) ** 0.5)
        print("corr between GCU length and the difference", corr[0,1])
        return
    def statistics(self,atgc):
        self.functions_names = ["GCU", "CGA"]
        self.l_functions = len(self.functions_names)
#        for atgc in self.ATGCs:  # for each atgc
        self.real_data(atgc)  # calculate statistics for the atgc
        self.print_res()
        return
get_args()
#ind = IndelRealDataStatistics()
