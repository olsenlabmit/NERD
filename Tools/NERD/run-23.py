from xlwt import Workbook 

from M23.Module_2_3 import format_bigsmiles, sub_obj_with_Bk, substructure_search_rxns, extract_arms, coarse_grain_arm
from polymersearch.search_tools import generate_NetworkX_graphs

f = open("Info-Drawing.txt", "r")
num_precursors = 0
inputs = []
for x in f:
    x = x.strip()
    if "NEW" in x:
        num_precursors += 1
    inputs.append(x) 

# The user can also encode their own information, 
# writing in the BigSMILES, all stochastic objects, 
# and chain length. For small molecules, the dictionary 
# will only have concentration and number {}
information = dict()
n = 0
for i in range(num_precursors):
    n += 2
    input_T = int(inputs[n])
    n += 2
    input_kmc = int(inputs[n])
    n += 2
    input_conversion = float(inputs[n])
    n += 2
    concentration = float(inputs[n])
    n += 2
    number_in_simulation = int(inputs[n])
    n += 2
    final_bigsmiles = inputs[n]
    information[final_bigsmiles] = dict()
    information[final_bigsmiles]["concentration"] = float(concentration)
    information[final_bigsmiles]["number"] = float(number_in_simulation)
    n += 2
    num_objects = int(inputs[n])
    n += 2
    if num_objects != 0:
        for j in range(num_objects):
            information[final_bigsmiles][j] = dict()
            information[final_bigsmiles][j]["object"] = inputs[n]
            n += 2
            information[final_bigsmiles][j]["Da"] = float(inputs[n])
            n += 2
            information[final_bigsmiles][j]["kDa"] = float(inputs[n])
            n += 2
            information[final_bigsmiles][j]["DP"] = float(inputs[n])
            n += 2
            information[final_bigsmiles][j]["Mo"] = float(inputs[n])
            n += 2
            information[final_bigsmiles][j]["b"] = float(inputs[n])
            n += 2
    num_precursors += 1

print(information)

# information = {"CC(C)(CCCCO)C(=O)O{[>][<]CCO[>][<]}C(=O)C(C)(C)CCCCO": 
#             {0: {'object': '{[][<]CCO[>][]}', 'Da': 0, 'kDa': 10, 'DP': 0, 'Mo': 137, 'b': 11}, 
#             "concentration": 20, "number": 20000},
#             'OCOCC(CO{[>][<]CC(C)[>][<]}CO)(COCO)C{[>][<]CCO[>][<]}OCC(=O)O': {
#             0: {'object': '{[][<]CC(C)[>][]}', 'Da': 0, 'kDa': 0, 'DP': 65, 'Mo': 137, 'b': 25},
#             1: {'object': '{[][<]CCO[>][]}', 'Da': 0, 'kDa': 0, 'DP': 265, 'Mo': 48, 'b': 8},
#             "concentration": 20, "number": 20000}
#             }

smiles = []
for bigsmiles in information:
    sm = format_bigsmiles(bigsmiles)
    
    # this makes it easier to visualize, since [Bk] is usually associated with stochastic objects
    sm = sub_obj_with_Bk(sm)[0]
    smiles.append(sm)

    # replace [Bk] with [Lr] since the graph generation code already uses [Bk]
    information[bigsmiles]["smiles"] = sm.replace("[Bk]", "[Lr]")

# reaction detection
matches = substructure_search_rxns(smiles) 

# Prompts the user to choose the reaction template for coarse graining
template_num = input('What reaction template do you choose from Info-Rxn?\n') 
matches = matches[int(template_num) - 1]

# extract and coarse grain each arm for N and b, determine largest b
count_bigsmiles = 0
b_max = 0
for bigsmiles in information:

    smiles_graphs = generate_NetworkX_graphs(information[bigsmiles]["smiles"]) 
    
    # extracts arms with stochastic objects exactly as they appear in the string
    # [arm smiles, A or B, stochastic object, is_linear]
    arms = extract_arms(bigsmiles, information[bigsmiles]["smiles"], smiles_graphs, matches[count_bigsmiles]) 
    count_bigsmiles += 1

    # first get all (N, b) pairs for all polymeric arms and determine the largest b
    count_object = 0
    pairs = []
    for attributes in arms:
        # if polymer
        if "Lr" in attributes[0]: 
            N_arm, b_arm = coarse_grain_arm(smiles=attributes[0], object=attributes[2], input_T=input_T, 
                                            information=information[bigsmiles][count_object], sub=True)
            pairs.append([N_arm, b_arm])
            count_object += 1
            if b_arm > b_max:
                b_max = b_arm

    information[bigsmiles]["arms"] = arms
    information[bigsmiles]["pairs"] = pairs

for bigsmiles in information:

    # stores topological distance from each A and B group to center of molecule
    arm_distances = [[],[]] 
    arms = information[bigsmiles]["arms"]
    pairs = information[bigsmiles]["pairs"]
    count_object = 0
    for attributes in arms:
        # if polymer
        if "Lr" in attributes[0]: 
            N_arm = pairs[count_object][0]
            b_arm = pairs[count_object][1]
            N_new = N_arm * b_arm**2 / b_max**2
            count_object += 1
            # if linker
            if attributes[3]:
                for i in range(2):
                    arm_distances[attributes[1]].append(N_new / 2)
            else:
                arm_distances[attributes[1]].append(N_new)
        else:
            if attributes[3]: 
                for i in range(2):
                    arm_distances[attributes[1]].append(0)
            else:
                arm_distances[attributes[1]].append(0)
    
    information[bigsmiles]["arm distance"] = arm_distances

file = open("Info-CG.txt", "w") 
file.write("Number of KMC" + "\n") 
file.write(str(input_kmc) + "\n") 
file.write("Conversion" + "\n") 
file.write(str(input_conversion) + "\n") 
file.write("Kuhn Length" + "\n")   
file.write(str(b_max) + "\n")
file.write("Number of Precursors" + "\n") 
file.write(str(len(information)) + "\n") 
for key in information:
    file.write("\nNEW" + "\n") 
    file.write("Concentration" + "\n") 
    file.write(str(information[key]["concentration"]) + "\n")
    file.write("Number in Simulation" + "\n") 
    file.write(str(information[key]["number"]) + "\n")
    N_template = information[key]["arm distance"]
    file.write("A" + "\n") 
    file.write(str(len(N_template[0])) + "\n")   
    for A in N_template[0]:
        file.write(str(A) + "\n") 
    file.write("B" + "\n")     
    file.write(str(len(N_template[1])) + "\n")
    for B in N_template[1]:
        file.write(str(B) + "\n") 
file.close()