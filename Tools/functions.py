import pandas as pd
import visualization
import sys, os

import polymersearch.search_tools
import polymersearch.fast 
import polymersimilarity.similarity_tools
from polymerscribe.generate import generate_data

# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore
def enablePrint():
    sys.stdout = sys.__stdout__

bigsmiles = "[H]{[>][<]C(C[>])c1ccccc1[<]}C(C)CC"

# Visualize BigSMARTS and BigSMILES
# graphs = polymersearch.search_tools.generate_NetworkX_graphs(bigsmiles, is_bigsmarts = False)
# graphs = visualization.set_colors(graphs, "atomistic")
# visualization.visualize_NetworkX_graphs(graphs, "atomistic", "atomistic_original") 
# graphs = visualization.set_colors(graphs, "multidigraph")
# visualization.visualize_NetworkX_graphs(graphs, "multidigraph", "multidigraph_original") 
# quit()

# Generate filters for BigSMARTS/BigSMILES
filters = polymersearch.fast.generate_filters(bigsmiles)
print("filters: ", filters)

# Validation
# filter_unit_testing()

# Fast Search
# bigsmarts = "{[][<]CCO[>][<]}CCO{[>][<]CCO[>][<]}CCO{[>][<]CCO[>][]}"
# bigsmiles = "{[][<]CCO[>][<]}CCO{[>][<]CCO[>][<]}CCO{[>][<]CCO[>][]}"
# query_graphs = polymersearch.search_tools.generate_NetworkX_graphs(bigsmarts, is_bigsmarts = True) 
# target_graphs = polymersearch.search_tools.generate_NetworkX_graphs(bigsmiles, is_bigsmarts = False) 
# x = polymersearch.fast.search_small_molecules(query_graphs, target_graphs) 

# generate polymerscribe data
# generate_data(to_excel=True)

# similarity
query = "{[][<]CCO[>][]}"
targets = [
            "{[][<]CCCO[>][]}",
            "{[][$]CC[$][]}",
            "{[][$]CC=CC[$][]}",
            "{[][<]NCC(=O)[>][]}",
            "{[][<]C(C)C(=O)O[>][]}",
            "{[][$]CC(C)(C(=O)OC)[$][]}",
            "{[][$]CC(OC(=O)C)[$][]}",
            "{[][$]CC(c1ccccc1)[$][]}",
            "{[][$]CC(C)(c1ccccc1)[$][]}",
            "{[][<]CCO[>],[<]CC(C)O[>][]}",
            "{[][$]CC(c1ccccc1)[$][$]}{[>][<]CCO[>][]}",
            "{[][$]CC(C(=O)O)[$],[$]CC(C(=O)O{[>][<]CCO[>][<]}C)[$][]}",
            "C([#Arm])([#Arm])([#Arm])[#Arm].{#Arm=CO{[<][>]CCO[<][>]}}",
            "{Cc1c([<1[<1]1])c([<1[<1]1])c(N(C)C)c3c1C4C2CC(C([>1[>1]2])C2[>1[>1]2])C34}",
            "C1C{[$][$]CC(c1ccccc1)[$][$]}CCCNC(=O)C1",
            "{[][<]NCCC{[>][<][Si](C)(C)O[>][<]}[Si](C)(C)CCCN[<],[>]C(=O)NC(CC1)CCC1CC(CC2)CCC2NC(=O)[>][]}",
            "{[][>]C(=O)Nc1ccc(C)c(c1)NC(=O)[>],[<]OCC{[>][<]OCC[>][<]}O[<],[<]OCCCO[<][]}",
            "{[][<]CCCCC(C)(C)C(=O)O{[>][<]CCO[>][<]}C(=O)C(C)(C)CCCC[<],[>]n1cc([<2])nn1,[>2]COCC(COC[>2])(COC[>2])COC[>2][]}"
        ] 

# blockPrint()
graphs = polymersearch.search_tools.generate_NetworkX_graphs(query, is_bigsmarts = False) 
v1 = polymersimilarity.similarity_tools.ECFP(graphs = graphs)
# enablePrint()
print(v1) 

for t in targets:
    graphs = polymersearch.search_tools.generate_NetworkX_graphs(t, is_bigsmarts = False) 
    v2 = polymersimilarity.similarity_tools.ECFP(graphs = graphs) 
    cosine, dice, manhattan, tanimoto = polymersimilarity.similarity_tools.pairwise_similarity(v1, v2) 
    enablePrint()
    print(cosine) 
    blockPrint()
quit()

# BigSMARTS Graph Traversal
# bigsmarts = "{[][]}"
# bigsmiles = "{[][<]OCCOC(=O)CCC(=O)[>][]}"
# found = polymersearch.search_tools.graph_traversal_match(bigsmarts, bigsmiles) 
# print(found) 

# Validation
# graph_traversal_unit_testing(bigsmarts_testing = [0,-1], sheet_name = "Extras") 