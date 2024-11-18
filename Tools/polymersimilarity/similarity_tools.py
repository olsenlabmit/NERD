from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
import numpy as np
import itertools
import networkx as nx    
import copy

def in_ring(atomistic, atom):
    symbols = nx.get_node_attributes(atomistic, "symbol")
    for neighbor in atomistic[atom]:
        count = 0
        for path in nx.all_simple_paths(atomistic, atom, neighbor):
            found = False
            for p in path:
                if ">" in symbols[p] or "<" in symbols[p] or "$" in symbols[p]:
                    found = True
                    break
            if not found:
                count += 1
            if count == 2:
                return True
    return False

def is_descriptor(symbol):
    return ">" in symbol or "<" in symbol or "$" in symbol

def bonds_compat_desc(list_bonds, bonds):
    path = []
    for i in range(len(list_bonds) - 1):
        bond_id = tuple(sorted([list_bonds[i], list_bonds[i + 1]]))
        path.append(bonds[bond_id][:2])

    possible = [
            ["1_", "2_"], 
            ["2_", "1_"], 
            ["1_", "2_", "1_", "2_"], 
            ["1_", "2_", "2_", "1_"], 
            ["2_", "1_", "1_", "2_"],
            ["2_", "1_", "2_", "1_"],
            ["1a", "1b"],
            ["1b", "1a"],  
            ["2a", "2b"], 
            ["2b", "2a"]
        ]
    
    return path in possible

def ECFP(graphs, num_iter = 2, fp_size = 1024):
    # https://chemicbook.com/2021/03/25/a-beginners-guide-for-understanding-extended-connectivity-fingerprints.html
    # remove [H] from graph

    atomistic = graphs["atomistic"]
    symbols = nx.get_node_attributes(atomistic, "symbol")
    for key in symbols:
        if symbols[key] == "H":
            atomistic.remove_node(key)
    symbols = nx.get_node_attributes(atomistic, "symbol")
    formal_charge = nx.get_node_attributes(atomistic, "formal_charge")
    valence = nx.get_node_attributes(atomistic, "valence")
    atomic_num = nx.get_node_attributes(atomistic, "atomic_num")
    mass = nx.get_node_attributes(atomistic, "mass")
    num_hs = nx.get_node_attributes(atomistic, "num_hs") 
    bonds = nx.get_edge_attributes(atomistic, "bond_type")

    all_identifiers = dict()
    feature_list = [] 
    for node in symbols:
        if not is_descriptor(symbols[node]) and symbols[node] != "":
            neighbors = len(atomistic[node]) 
            valence_atoms = valence[node]
            hydrogens_atom = num_hs[node]
            atomic_atom = atomic_num[node]
            mass_atom = mass[node]
            charge_atom = formal_charge[node] 
            ring = in_ring(atomistic, node)
            
            print((neighbors, valence_atoms - hydrogens_atom, atomic_atom, mass_atom, charge_atom, hydrogens_atom, ring))
            identifier = hash((neighbors, valence_atoms - hydrogens_atom, atomic_atom, mass_atom, charge_atom, hydrogens_atom, ring))
            print(identifier)
            all_identifiers[(node, 0)] = [identifier]
            for identifier in all_identifiers[(node, 0)]:
                feature_list.append(identifier)
    
    iteration = 1
    while iteration <= num_iter:
        print("Iteration: ", iteration)
        for node in symbols:
            if is_descriptor(symbols[node]) or symbols[node] == "":
                continue   
            
            print("node: ", node)
            new_identifiers = []
            for atom_hash in all_identifiers[(node, iteration - 1)]:
                original = (iteration, atom_hash)
                print("original: ", original)

                neighbor_hash = []
                for neighbor_1 in atomistic[node]:
                    if is_descriptor(symbols[neighbor_1]):
                        continue

                    n = [] 
                    for bond in bonds:
                        if bond == tuple(sorted([node, neighbor_1])):
                            if "SINGLE" in bonds[bond]:
                                order = 1
                            elif "DOUBLE" in bonds[bond]:
                                order = 2
                            elif "TRIPLE" in bonds[bond]:
                                order = 3 
                            else:
                                order = 4
                    for n_hash in all_identifiers[(neighbor_1, iteration - 1)]:
                        n.append((order, n_hash))
                    neighbor_hash.append(n) 
                print("neighbor hash: ", neighbor_hash)
                    
                across_hash = [] 
                for neighbor_1 in atomistic[node]:
                    if not is_descriptor(symbols[neighbor_1]):
                        continue
                    for neighbor_2 in atomistic[neighbor_1]:
                        if symbols[neighbor_2] == "":
                            for neighbor_3 in atomistic[neighbor_2]:
                                for neighbor_4 in atomistic[neighbor_3]:
                                    for bond in bonds:
                                        if bond == tuple(sorted([neighbor_3, neighbor_4])):
                                            if "SINGLE" in bonds[bond]:
                                                order = 1
                                            elif "DOUBLE" in bonds[bond]:
                                                order = 2
                                            elif "TRIPLE" in bonds[bond]:
                                                order = 3 
                                            else:
                                                order = 4
                                    list_bonds = [node, neighbor_1, neighbor_2, neighbor_3, neighbor_4]
                                    a = symbols[node] != "" and not is_descriptor(symbols[node])
                                    b = is_descriptor(symbols[neighbor_1])
                                    c = symbols[neighbor_2] == ""
                                    d = is_descriptor(symbols[neighbor_3])
                                    e = symbols[neighbor_4] != "" and not is_descriptor(symbols[neighbor_4])
                                    
                                    compatible = bonds_compat_desc(list_bonds, bonds)
                                    repeat_nodes = len(set(list_bonds)) == 5
                                    if not(a and b and c and d and e and compatible and repeat_nodes): 
                                        continue
                                    for a_hash in all_identifiers[(neighbor_4, iteration - 1)]:
                                        across_hash.append((order, a_hash))
                        else:
                            for bond in bonds:
                                if bond == tuple(sorted([neighbor_1, neighbor_2])):
                                    if "SINGLE" in bonds[bond]:
                                        order = 1
                                    elif "DOUBLE" in bonds[bond]:
                                        order = 2
                                    elif "TRIPLE" in bonds[bond]:
                                        order = 3 
                                    else:
                                        order = 4
                            list_bonds = [node, neighbor_1, neighbor_2]
                            if not bonds_compat_desc(list_bonds, bonds):
                                continue
                            for a_hash in all_identifiers[(neighbor_2, iteration - 1)]:
                                across_hash.append((order, a_hash))
                print("across hash: ", across_hash)

                list_length = []
                for n in neighbor_hash:
                    list_length.append(list(range(len(n))))
                combos = [p for p in itertools.product(*list_length)] 

                if across_hash:
                    for c in combos: 
                        for a in across_hash:
                            neighbors = []
                            for i in range(len(neighbor_hash)):
                                neighbors.append(neighbor_hash[i][c[i]])
                            print("neighbors: ", neighbors)
                            neighbors.append(a) 
                            print("neighbors + across: ", neighbors)
                            neighbors = sorted(neighbors)
                            print("sorted:", neighbors)

                            proper_list = []
                            for item in original:
                                proper_list.append(item)
                            for pair in neighbors:
                                for item in pair:
                                    proper_list.append(item)
                            print("proper list: ", proper_list)
                            identifier = hash(tuple(proper_list)) 
                            print("identifier: ", identifier)

                            new_identifiers.append(identifier) 
                            feature_list.append(identifier) 
                else:
                    for c in combos: 
                        neighbors = []
                        for i in range(len(neighbor_hash)):
                            neighbors.append(neighbor_hash[i][c[i]])
                        print("neighbors: ", neighbors)                       
                        neighbors = sorted(neighbors)
                        print("sorted neighbors: ", neighbors)

                        proper_list = []
                        for item in original:
                            proper_list.append(item)
                        for pair in neighbors:
                            for item in pair:
                                proper_list.append(item)
                        print("proper list: ", proper_list)
                        identifier = hash(tuple(proper_list)) 
                        print("identifier: ", identifier)

                        new_identifiers.append(identifier) 
                        feature_list.append(identifier) 
            
            all_identifiers[(node, iteration)] = new_identifiers
            print("Updated: ", node, iteration, new_identifiers)
    
        iteration += 1
        # quit()

    print(feature_list)

    # convert to an array
    fp = np.zeros(fp_size)
    remainders = []
    for f in feature_list:
        remainders.append(f % fp_size)
    print(remainders)
    for r in remainders:
        fp[r] = 1
    print(fp)
    
    return fp    

def atom_pair(graphs, fp_size = 1024):
    # https://pubs.acs.org/doi/epdf/10.1021/ci00046a002
    
    atomistic = graphs["atomistic"]
    symbols = nx.get_node_attributes(atomistic, "symbol")
    is_aromatic = nx.get_node_attributes(atomistic, "is_aromatic")
    atomic_num = nx.get_node_attributes(atomistic, "atomic_num")

    feature_list = []
    explored = []
    for node_1 in symbols:
        for node_2 in symbols:
            a = not (is_descriptor(symbols[node_1]) or symbols[node_1] == "")
            b = not (is_descriptor(symbols[node_2]) or symbols[node_2] == "")
            c = node_1 != node_2
            d = sorted([node_1, node_2]) not in explored
            if a and b and c and d:
                explored.append(sorted([node_1, node_2]))
                path = nx.shortest_path(atomistic, node_1, node_2)

                def atom_type(node):
                    atomic_atom = atomic_num[node]
                    is_a = is_aromatic[node]
                    return (atomic_atom, is_a)
                
                f1 = atom_type(node_1)
                f2 = atom_type(node_2)
                f = sorted([f1, f2])
                b = tuple(f + [len(path)])
                feature_list.append(hash(b))

    fp = np.zeros(fp_size)
    remainders = []
    for f in feature_list:
        remainders.append(f % fp_size)
    for r in remainders:
        fp[r] = 1

    return fp

def topological_torsion(graphs, fp_size = 1024):
    atomistic = graphs["atomistic"]
    symbols = nx.get_node_attributes(atomistic, "symbol")
    bonds = nx.get_edge_attributes(atomistic, "bond_type")
    is_aromatic = nx.get_node_attributes(atomistic, "is_aromatic")
    atomic_num = nx.get_node_attributes(atomistic, "atomic_num")

    all_paths = []
    def dfs(current, parent, path):
        if len(path) == 4:
            all_paths.append(copy.deepcopy(path))
            path.pop()
            return
        neighbors = atomistic[current]
        for neighbor_1 in neighbors:
            if is_descriptor(symbols[neighbor_1]):
                for neighbor_2 in atomistic[neighbor_1]:
                    list_bonds = [current, neighbor_1, neighbor_2]
                    if bonds_compat_desc(list_bonds, bonds):
                        if neighbor_2 != current and neighbor_2 != parent:
                            path.append(neighbor_2)
                            dfs(neighbor_2, current, path)
            else:
                if neighbor_1 != parent:
                    path.append(neighbor_1)
                    dfs(neighbor_1, current, path)

    for key in symbols:
        if not (is_descriptor(symbols[key]) or symbols[key] == ""):
            dfs(key, -1, [key]) 
    print(all_paths) 

    feature_list = []
    for path in all_paths:
        to_hash = []
        def atom_type(node):
            atomic_atom = atomic_num[node]
            is_a = is_aromatic[node]
            return (atomic_atom, is_a)
        for p in path:
            to_hash.append(atom_type(p))
        feature_list.append(hash(tuple(to_hash)))   

    fp = np.zeros(fp_size)
    remainders = []
    for f in feature_list:
        remainders.append(f % fp_size)
    for r in remainders:
        fp[r] = 1

    return fp

def pairwise_similarity(x, y):
    # https://docs.eyesopen.com/toolkits/python/graphsimtk/measure.html

    A = 0
    B = 0
    onlyA = 0
    onlyB = 0
    bothAB = 0
    neitherAB = 0
    for key in range(len(x)):
        if x[key] == 1:
            A += 1
        if y[key] == 1:
            B += 1
        if x[key] == 1 and y[key] == 0:
            onlyA += 1
        if x[key] == 0 and y[key] == 1:
            onlyB += 1
        if x[key] == 1 and y[key] == 1:
            bothAB += 1
        if x[key] == 0 and y[key] == 0:
            neitherAB += 1
    
    cosine = bothAB/np.sqrt(A*B)
    dice = 2 * bothAB / (A + B)
    manhattan = (onlyA + onlyB) / len(x)
    tanimoto = bothAB / (onlyA + onlyB + bothAB)
    
    return cosine, dice, manhattan, tanimoto