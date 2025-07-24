import argparse

from scipy.optimize import fsolve
def func(x,p):
    theta = x[0] #660
    alph =  x[1] #0.675
    n_genes = p[0][0]
    cogs_final = p[1][0]
    cogs_of_one_final = p[2][0]
    cogs =1
    cogs_of_one = 1

    for i in range(n_genes):
        add_tn = (theta + alph * cogs) / (theta + i)
        cogs += add_tn
        cogs_of_one += add_tn - (1 - alph) * cogs_of_one / (theta + i)
    return [cogs - cogs_final , cogs_of_one - cogs_of_one_final]

def compute_parameters(p):
    root = fsolve(func,[2000,0.5] ,p)
    return root
#print('Expected number of cogs:', CogE)


parser = argparse.ArgumentParser(description ='compute the parameters for chinese restaurant.')
parser.add_argument('genes',
                    metavar ='genes',
                    type = int,
                    nargs = 1,
                    help ='number of genes in the genome')

parser.add_argument('cogs',
                    metavar ='cogs',
                    type = float,
                    nargs = 1,
                    help ='expected number of cogs')

parser.add_argument('one_gene_cogs',
                    metavar ='ogc',
                    type = float,
                    nargs = 1,
                    help ='expected number of one-gene cogs')

parser.add_argument(dest ='cp',
                    action ='store_const',
                    const = compute_parameters,
                    help ='compute the parameters theta and alpha')

args = parser.parse_args()
print(args.cp([args.genes,args.cogs,args.one_gene_cogs]))