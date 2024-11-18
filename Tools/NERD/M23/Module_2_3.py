import rdkit
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem, Descriptors, rdChemReactions, Draw, AllChem
import pandas as pd
import re
import copy
import networkx as nx
import chemprop

from polymersearch.graphs import get_objects, get_repeats, desc_regex
from polymersearch.search_tools import graph_traversal_match

def sub_obj_with_Bk(bigsmiles):
    smiles = ""
    counter = 0
    objects = get_objects(bigsmiles)[0]
    indices = get_objects(bigsmiles)[1]
    for i in range(len(objects)):
        smiles += bigsmiles[counter:indices[i]]
        smiles += "[Bk]"
        counter += (indices[i] - counter) + len(objects[i])
    smiles += bigsmiles[counter:]
    return smiles, objects

def contains_both(precursors, A, B):
    containsA = False
    containsB = False
    for precursor in precursors:
        if precursor.HasSubstructMatch(A):
            containsA = True
        if precursor.HasSubstructMatch(B):
            containsB = True
    return containsA and containsB

def substructure_search_rxns(smiles):
    precursors = []
    for i in range(len(smiles)):
        precursors.append(Chem.MolFromSmiles(smiles[i]))

    all_atoms = []
    rxn_split = []
    databases = ["JCIM_2011", "NERD"]
    for d in databases:
        df = pd.read_excel('M23/Database.xlsx', sheet_name = d) 
        SMARTS = df['SMARTS'] 
        for rxn in SMARTS:
            s = rxn.split(">")
            reactants = s[0].split(".")
            for i in range(len(reactants)):
                for j in range(i + 1, len(reactants)):
                    A = Chem.MolFromSmarts(reactants[i])
                    B = Chem.MolFromSmarts(reactants[j])
                    r_split = []
                    r_all = []
                    if contains_both(precursors, A, B):
                        for molecule in precursors:
                            r_split.append([0,0])
                            matches = molecule.GetSubstructMatches(A)
                            h = [element for tupl in matches for element in tupl]
                            r_split[-1][0] = sorted(matches)
                            matches = molecule.GetSubstructMatches(B)
                            h += [element for tupl in matches for element in tupl]
                            r_all.append(h)
                            r_split[-1][1] = sorted(matches)
                    if r_split:
                        rxn_split.append(r_split)
                        all_atoms.append(r_all)

    if rxn_split == []:
        return "No Reaction!"
    
    molecules = []
    h = []
    legend = []
    for i in range(len(all_atoms)):
        for p in range(len(precursors)):
            molecules.append(precursors[p])
            h.append(all_atoms[i][p])
            legend.append('Template #' + str(i + 1))

    if all_atoms:
        img = Draw.MolsToGridImage(molecules, legends = legend, subImgSize = (700, 700), molsPerRow = len(precursors) * 2, highlightAtomLists = h)
        img.save('Info-Rxn.png')
    
    return rxn_split

def format_bigsmiles(bigsmiles):
    replace = {";!{[][]}":"","!{[][]}":"","[H]{":"{","}[H]":"}","[<]":"[<1]","[>]":"[>1]","[$]":"[$1]","?*":"[Fm][Md]","[<=]":"[<=1]","[>=]":"[>=1]","[$=]":"[$=1]"} 
    for key in replace:
        bigsmiles = bigsmiles.replace(key, replace[key])
    return bigsmiles

def bigsmiles_to_repeat(bigsmiles):
    objects = get_objects(bigsmiles)[0]
    for o in objects:
        repeats = get_repeats(o)[0]
        for r in repeats:
            repeat = r
    return repeat

def replace_descriptors_Cf(repeat):
    descriptor_locations = [(d.start(0), d.end(0)) for d in re.finditer(desc_regex, repeat)] 
    descriptors = []
    for d in descriptor_locations:
        descriptors.append(repeat[d[0]:d[1]])
    for d in descriptors:
        repeat = repeat.replace(d,"[Cf]")
    return repeat

def extract_arms(bigsmiles, smiles, graphs, matches):

    # figure out bonds to break, if there is exactly one path between atoms, then break
    atomistic = graphs["atomistic"]
    atomistic_symbols = nx.get_node_attributes(atomistic, "symbol")
    atomistic_bonds = nx.get_edge_attributes(atomistic, "bond_type")
    break_bond = []
    for key in atomistic_bonds:
        path_found = []
        for path in nx.all_simple_paths(atomistic, key[0], key[1]):
            desc_found = False
            for node in path:
                if "<" in atomistic_symbols[node] or ">" in atomistic_symbols[node] or "$" in atomistic_symbols[node]:
                    desc_found = True
                    break
            if not desc_found:                            
                path_found.append(path)
        if len(path_found) == 1:
            break_bond.append(key)    
    # print("break_bond: ", break_bond) 

    topology = graphs["topology"]
    topology_symbols = nx.get_node_attributes(topology, "symbol")
    topology_ids = nx.get_node_attributes(topology, "ids") 
    non_descriptors = []
    for key in topology_symbols:
        if not ("<" in topology_symbols[key] or ">" in topology_symbols[key] or "$" in topology_symbols[key]):
            non_descriptors.append(topology_ids[key])
    
    # get all paths between ends of molecule, this is the backbone
    start_ends = []
    for m in matches:
        for group in m:
            start_ends.append(group[0]) 
    shortest_paths = []
    for i in range(len(start_ends)):
        for j in range(i + 1, len(start_ends)):
            for path in nx.all_simple_paths(atomistic, start_ends[i], start_ends[j]):
                shortest_paths.append(path)
                break
    # print("Paths Between Descriptors: ", shortest_paths)

    # only break bonds along the backbone paths
    break_bond_modified = []
    for bond in break_bond:
        found = False
        for path in shortest_paths:
            if bond[0] in path and bond[1] in path:
                found = True
                break
        if found:
            break_bond_modified.append(bond)
    # print("Bonds to Break: ", break_bond_modified) 
    
    # add states along the shortest path
    state_machine = copy.deepcopy(atomistic)
    state = max(list(nx.get_node_attributes(state_machine, "symbol").keys()))
    for breaking in break_bond_modified:
        state += 1
        state_machine.remove_edge(breaking[0], breaking[1])
        state_machine.add_node(state)
        state_machine.add_edge(breaking[0], state) 
        state_machine.add_edge(state, breaking[1]) 
    graphs["atomistic"] = state_machine 
    
    # determine the atoms in the branch point alphabet
    branched_atoms = []
    for atom in atomistic_symbols:
        # get all atoms in the same alphabet, 
        # key and atom are in the same alphabet 
        # if there is a path with no state in between
        branch = [atom]
        for key in atomistic_symbols:
            if "<" in atomistic_symbols[key] or ">" in atomistic_symbols[key] or "$" in atomistic_symbols[key]: 
                continue 
            for path in nx.all_simple_paths(state_machine, key, atom):
                valid = True
                for p in path:
                    if p not in atomistic_symbols:
                        valid = False
                        break
                if valid:
                    branch.append(key)
                    break
        # determine the number of neighboring states
        neigh_states = set()
        for g in branch:
            for n in state_machine[g]:
                if n not in atomistic_symbols:
                    neigh_states.add(n)
        # if this number is greater than two, 
        # then you have the branch point alphabet
        if len(neigh_states) > 2:
            branched_atoms = branch
            break

    # extract atoms in the same arm in the order they appear in the string
    # if there is a path that crosses the branch point,
    # then the atoms are not a part of the same arm
    arms = []
    for atom in atomistic_symbols:
        arm_groups = [atom]
        for key in atomistic_symbols:
            for path in nx.all_simple_paths(graphs["atomistic"], key, atom):
                valid = True
                for p in path:
                    if p in branched_atoms:
                        valid = False
                        break
                if valid:
                    arm_groups.append(key)
                    break
        if sorted(arm_groups) not in arms and sorted(arm_groups) != sorted(branched_atoms):
            arms.append(arm_groups)

     # this gets a list of the stochastic objects 
     # in the order they appear in the string
    repeat_units = sub_obj_with_Bk(bigsmiles)[1]

    # stores in each arm [arm smiles, A or B group, stochastic object or None, is_linear]
    if branched_atoms == []: # if linear
        for a in arms: 
            for A_or_B in range(len(matches)):
                for group in matches[A_or_B]: 
                    if group[0] in a:
                        # which set of matches can be found in that arm
                        if "[Lr]" in smiles:
                            # gets the first and only stochastic object
                            return [[smiles, A_or_B, repeat_units[0], True]] 
                        else:
                            return [[smiles, A_or_B, None, True]]
        
    else: # if branched 
        count_branch = 0
        information = [] 
        for a in arms: 
            found = False
            # returns (smiles, root atom)
            extracted = chemprop.interpret.extract_subgraph(smiles, a) 
            # gets a single atom from the extracted list
            if extracted[1][0] not in branched_atoms: 
                for A_or_B in range(len(matches)):
                    for group in matches[A_or_B]: 
                        if group[0] in a:
                            found = True
                            if "[Lr]" in extracted[0]:
                                information.append([extracted[0], A_or_B, repeat_units[count_branch], False]) 
                                count_branch += 1
                            else:
                                information.append([extracted[0], A_or_B, None, False]) 
                            break
                    if found:
                        break

    return information

def coarse_grain_arm(smiles, object, input_T, information, sub):

    stochastic_object = format_bigsmiles(object)
    repeat = bigsmiles_to_repeat(stochastic_object)
    input_repeat_Cf = replace_descriptors_Cf(repeat)
    num = Chem.MolFromSmiles(input_repeat_Cf).GetNumHeavyAtoms() - 2
    
    def get_Mo_b(): 
        data = pd.read_excel("M23/Database.xlsx", sheet_name = "Kuhn")
        BigSMILES = data["BigSMILES"]
        T = data["T (K)"]
        Mk = data["Mo (g/mol)"]
        k = data["b (A)"]
        found = []
        def sim():
            for i in range(len(BigSMILES)):
                db_r = format_bigsmiles(BigSMILES[i])
                db_r = bigsmiles_to_repeat(db_r)
                db_r = replace_descriptors_Cf(db_r)
                db_r = db_r.replace("[Cf]", "[*]")
                r2 = input_repeat_Cf.replace("[Cf]", "[*]")
                a = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(db_r), radius = 2, useChirality = True, nBits=2048)
                b = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(r2), useChirality = True, radius = 2, nBits=2048)
                s = DataStructs.TanimotoSimilarity(a, b)
                found.append((BigSMILES[i], (1 - s), abs(T[i] - input_T), Mk[i], k[i]))
        if sub:
            for i in range(len(BigSMILES)):
                if graph_traversal_match(stochastic_object, BigSMILES[i]):
                    db_r = format_bigsmiles(BigSMILES[i])
                    db_r = bigsmiles_to_repeat(db_r)
                    db_r = replace_descriptors_Cf(db_r)
                    num2 = Chem.MolFromSmiles(db_r).GetNumHeavyAtoms() - 2
                    found.append((BigSMILES[i], num2 - num, abs(T[i] - input_T), Mk[i], k[i]))
            if not found:
                sim()
        else:
            sim()
        
        found.sort(key=lambda x: x[1], reverse=False)
        f = [found[0]]
        i = 1
        target = found[0][1]
        while i < len(found) and found[i][1] == target:
            f.append(found[i])
            i += 1
        f.sort(key=lambda x: x[2], reverse=False)
        Mk = f[0][3]
        k = f[0][4]
        return Mk, k

    if information["Mo"] == 0:
        Mo, b_arm = get_Mo_b()
    else:
        # for now, the user enters the information
        Mo = information["Mo"]
        b_arm = information["b"]

    if information["Da"] != 0:
        molar_mass = information["Da"]
    elif information["kDa"] != 0:
        molar_mass = information["kDa"] * 1000
    else:
        # add endgroup minus Lr
        # add repeat unit * DP
        rdkit_smiles = Chem.MolFromSmiles(smiles)
        molar_mass = rdkit.Chem.Descriptors.ExactMolWt(rdkit_smiles) - rdkit.Chem.Descriptors.ExactMolWt(Chem.MolFromSmiles("[Lr]"))
        molar_mass += (rdkit.Chem.Descriptors.ExactMolWt(Chem.MolFromSmiles(input_repeat_Cf)) - 2 * rdkit.Chem.Descriptors.ExactMolWt(Chem.MolFromSmiles("[Cf]"))) * information["DP"] 
    return molar_mass / Mo, b_arm