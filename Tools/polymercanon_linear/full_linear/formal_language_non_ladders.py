import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import re
import rdkit
from rdkit import Chem
import pandas as pd
import pylab
import numpy as np
import priority_rules_non_ladders  
import chemprop
from pyvis.network import Network
from rdkit.Chem import Draw
import copy
import matplotlib.pyplot as plt
import pydot

from BigSMILES_parser.BigSMILES_BigSmilesObj import BigSMILES

def get_compatible_descriptor(d):
    if "<" in d: return d.replace("<", ">")
    elif ">" in d: return d.replace(">", "<")
    else: return d

def process(bigsmiles):
    bigsmiles = priority_rules_non_ladders.process(bigsmiles)
    d = []
    d += re.findall(r"\[.\d+\]", bigsmiles)
    if bigsmiles.find("{[]") >= 0:
        bigsmiles = bigsmiles.replace("{[]", "{" + get_compatible_descriptor(d[0]))
    if bigsmiles.find("[]}") >= 0:
        bigsmiles = bigsmiles.replace("[]}", get_compatible_descriptor(d[-1])+"}")
    return bigsmiles

def parse_objects(bigsmiles):
    repeat_units = []
    terminals = []
    objects = priority_rules_non_ladders.get_objects(bigsmiles)[0]
    for object in objects:
        repeat_in_object = []
        repeats = priority_rules_non_ladders.get_repeats(object)
        for repeat in repeats:
            nested = priority_rules_non_ladders.get_objects(repeat)[0] 
            if len(nested) > 0:
                repeat_in_object.append(repeat)
            else:
                repeat = str(BigSMILES(repeat)).replace("-", "")
                repeat_in_object.append(repeat)
        repeat_units.append(repeat_in_object)
        left_terminal = object[1:object.find("]")+1]
        right_terminal = object[object.rfind("["):-1]
        terminals.append([left_terminal, right_terminal])
    return repeat_units, terminals

def parse_endgroups(bigsmiles):
    def extract_atoms(tree, absolute_start, extracted, start_atom, stop):
        extracted.add(start_atom)
        next_atoms = []
        for i in tree:
            if i[0] == start_atom and i[1] > absolute_start:
                extracted.add(i[1])
                if i[1] != stop:
                    next_atoms.append(i[1])

        for i in next_atoms:
            extracted = extract_atoms(tree, absolute_start, extracted, i, stop)
        return extracted
    
    # Replace all objects with Bk
    objects = priority_rules_non_ladders.get_objects(bigsmiles)
    bs = ""
    counter = 0
    for i in range(len(objects[0])):
        bs += bigsmiles[counter:objects[1][i]] + "[C:2][Bk][C:1]"
        counter = counter + len(objects[0][i]) + (objects[1][i] - counter)
    bs += bigsmiles[counter:]
    bigsmiles_cg = "[Bk][C:1]" + bs + "[C:2]"
    # print("after replacing objects: ", bigsmiles_cg)
    
    # Get all atoms between the stochastic objects, labeled Bk
    index = []
    m = Chem.MolFromSmiles(bigsmiles_cg)
    for q, atom in enumerate(m.GetAtoms()):
        if atom.GetSymbol() == "Bk":
            index.append(atom.GetIdx())
    x = RDKit_to_networkx_graph(Chem.MolFromSmiles(bigsmiles_cg))
    index.append(1000)
    end_groups = []
    for i in range(len(index) - 1):
        T = nx.dfs_tree(x, source=index[i])
        extracted = extract_atoms(list(T.edges), index[i], set(), index[i], index[i+1])
        end_group = chemprop.interpret.extract_subgraph(bigsmiles_cg, extracted)[0]
        end_group = end_group.replace("[Bk:1]","")
        end_group = end_group.replace("[Bk]","")
        end_group = end_group.replace("[C:1]","[<]")
        end_group = end_group.replace("[C:2]","[>]")
        # print("eg before: ", end_group)
        end_group = str(BigSMILES(end_group))
        # print("eg after: ", end_group)
        if end_group.find("[>1]") == 0:
            end_group = str(BigSMILES(end_group).writeStandard(noBondDesc=True,forward = False))
        else:
            end_group = str(BigSMILES(end_group).writeStandard(noBondDesc=True,forward = True))
        end_groups.append(end_group)
    
    return end_groups

def RDKit_to_networkx_graph(mol):
    G = nx.Graph()
    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                    bond.GetEndAtomIdx(),
                    bond_type=bond.GetBondType())
    return G

def generate_alphabets(cycle):
    cycle = "[U]" + cycle + "[U]"

    objects = priority_rules_non_ladders.get_objects(cycle)
    bs = ""
    counter = 0
    for i in range(len(objects[0])):
        bs += cycle[counter:objects[1][i]] + "[Bk]"
        counter = counter + len(objects[0][i]) + (objects[1][i] - counter)
    bs += cycle[counter:]
    cycle = bs
    # print("cycle: ", cycle)

    mol = Chem.MolFromSmiles(cycle)
    G = RDKit_to_networkx_graph(mol)
    path = nx.shortest_path(G, source=0, target=G.number_of_nodes()-1)
    # print("path: ", path)
    # check if object is in path
    counter = 0
    for q, atom in enumerate(mol.GetAtoms()):
        if atom.GetSymbol() == "Bk" and counter in path:
            return True
        counter += 1

    def sort_matches(a):
        x = []
        for i in range(len(a)):
            x.append(sorted(list(a[i])))
        return x
    
    break_at = []

    aliphatic_aliphatic = Chem.MolFromSmarts("A-A")
    aliphatic_aliphatic = mol.GetSubstructMatches(aliphatic_aliphatic)
    aliphatic_aliphatic = sort_matches(aliphatic_aliphatic)
    
    ringring = Chem.MolFromSmarts("[*;R][*;R]")
    ringring = sorted(mol.GetSubstructMatches(ringring))
    ringring = sort_matches(ringring)

    nonring_ring = Chem.MolFromSmarts("[*;!R][*;R]")
    nonring_ring = sorted(mol.GetSubstructMatches(nonring_ring))
    nonring_ring = sort_matches(nonring_ring)

    for t in aliphatic_aliphatic:
        if t not in ringring:
            break_at.append(t)

    for t in nonring_ring:
        break_at.append(t)
    
    for t in ringring:
        x = nx.all_simple_paths(G, source=t[0], target=t[1])
        if len(list(x)) == 1:
            break_at.append(t)
    #print("break_at", break_at)

    final = []
    for p in break_at:
        if p[0] in path and p[1] in path:
            final.append(sorted(p))   
    final.sort(key=lambda y: y[0])
    #print("final", final)

    def __extract_subgraph(mol, selected_atoms, dir1, dir2):
        """
        Extracts a subgraph from an RDKit molecule given a set of atom indices.

        :param mol: An RDKit molecule from which to extract a subgraph.
        :param selected_atoms: The atoms which form the subgraph to be extracted.
        :return: A tuple containing an RDKit molecule representing the subgraph
                 and a list of root atom indices from the selected indices.
        """
        selected_atoms = set(selected_atoms)
        roots = []
        for idx in selected_atoms:
            atom = mol.GetAtomWithIdx(idx)
            bad_neis = [y for y in atom.GetNeighbors() if y.GetIdx() not in selected_atoms]
            if len(bad_neis) > 0:
                roots.append(idx)

        new_mol = Chem.RWMol(mol)
        
        # Nathan
        if len(roots) == 1:
            roots.append(roots[0])
        left = Chem.rdchem.Atom("*")
        left.SetAtomMapNum(1)
        right = Chem.rdchem.Atom("*")
        right.SetAtomMapNum(2)
        new_mol.ReplaceAtom(dir1, left)
        new_mol.ReplaceAtom(dir2, right)
        selected_atoms.add(dir1)
        selected_atoms.add(dir2)

        for atom_idx in roots:
            atom = new_mol.GetAtomWithIdx(atom_idx)
            #atom.SetAtomMapNum(1)
            aroma_bonds = [bond for bond in atom.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.AROMATIC]
            aroma_bonds = [bond for bond in aroma_bonds if
                           bond.GetBeginAtom().GetIdx() in selected_atoms and bond.GetEndAtom().GetIdx() in selected_atoms]
            if len(aroma_bonds) == 0:
                atom.SetIsAromatic(False)

        remove_atoms = [atom.GetIdx() for atom in new_mol.GetAtoms() if atom.GetIdx() not in selected_atoms]
        remove_atoms = sorted(remove_atoms, reverse=True)
        for atom in remove_atoms:
            new_mol.RemoveAtom(atom)
        
        x = new_mol.GetMol()
        return Chem.MolToSmiles(x), roots

    series = []
    for f in range(len(final) - 1):
        start = final[f][1] # Get all atoms between 2nd atom of the previous pair and the 2nd atom of the next pair
        stop = final[f+1][1] # 2nd atom of the next pair
        s = list(range(start, stop))
        dir1 = final[f][0]
        dir2 = final[f+1][1]
        backbone = __extract_subgraph(Chem.MolFromSmiles(cycle),s,dir1,dir2)[0]
        mol = Chem.MolToSmiles(Chem.MolFromSmiles(backbone))
        if mol == "[Cf]([*:1])[*:2]": # ...except the empty string is between objects
            series.append("[*:1][*:2]")
        elif mol != "[*:1][*:2]" and mol != "[*:2][*:1]": # Remove empty strings
            series.append(mol)
    return series

def reverse_alphabets(bb):
    k = copy.deepcopy(bb)
    k.reverse()
    for i in range(len(k)):
        k[i] = k[i].replace("[*:1]","[Cf]")
        k[i] = k[i].replace("[*:2]","[*:1]")
        k[i] = k[i].replace("[Cf]","[*:2]")
    return k

def isomers(repeat):
    counter = 0
    alt = []
    while counter < len(repeat):
        string = ""
        if repeat[counter] == "C":
            string += repeat[counter]
            r = counter + 1
            parenthesis = 0
            if r < len(repeat) and repeat[r] == "(":
                parenthesis += 1
                string += repeat[r]
            r += 1
            while parenthesis != 0 and r < len(repeat):
                string += repeat[r]
                if repeat[r] == "(":
                    parenthesis += 1
                elif repeat[r] == ")":
                    parenthesis -= 1
                r += 1
            alt.append(string)
        counter += 1
    for i in range(len(alt) - 1):
        x = alt[i] + "=" + alt[i+1]
        ind = repeat.find(x)
        if ind != -1 and repeat[ind-1] not in ["\\","/"]:
            cis = "\\" + alt[i] + "=" + alt[i+1] + "/"
            trans = "/" + alt[i] + "=" + alt[i+1] + "/"
            repeat_cis = repeat[0:ind] + cis + repeat[ind+len(x):]
            repeat_trans = repeat[0:ind] + trans + repeat[ind+len(x):]
            try:
                img = Draw.MolsToGridImage([Chem.MolFromSmiles(repeat_cis)], subImgSize=(200, 200))
                img = Draw.MolsToGridImage([Chem.MolFromSmiles(repeat_trans)], subImgSize=(200, 200))
                img.save("static/test.jpeg")
                return [repeat_cis, repeat_trans]
            except:
                continue
    return [repeat]

def generate_paths(bigsmiles, left_desc, right_desc, fragments, backbone, object_id_highest):   

    print("bigsmiles enter function: ", bigsmiles, left_desc, right_desc)

    repeats, terminals = parse_objects(bigsmiles)
    end_group = parse_endgroups(bigsmiles)    

    for i in range(len(end_group)):
        if end_group[i] == "":
            end_group[i] = "[Cf]"

    origin = str(object_id_highest - 1)
    current = origin
    to_automata = []
    for i in range(len(end_group)):  

        prev = current
        current = str(object_id_highest)

        if i < len(end_group) - 1:
            
            # treat endgroups
            if i == 0:
                to_automata.append([left_desc + origin, end_group[i], terminals[i][0] + current])

            else:
                to_automata.append([terminals[i - 1][1] + prev, end_group[i], terminals[i][0] + current]) 

            backbone[end_group[i]] = generate_alphabets(end_group[i])
            backbone[end_group[i] + "_reverse"] = reverse_alphabets(backbone[end_group[i]])

            # treat repeat units
            for r in repeats[i]:
                smiles = r[4:-4]
                left = r[:4]
                right = r[-4:]
                wct = isomers(smiles)
                
                for k in wct:

                    print("repeat unit: ", k) 
                    object_in_backbone = generate_alphabets(k)

                    # if repeat unit has an object in the backbone
                    if object_in_backbone == True:

                        fragments, backbone, object_id_highest = generate_paths(k, left, right, fragments, backbone, int(current) + 1) 
                    
                    else:

                        to_automata.append([left + current, k, right + current])
                        backbone[k] = generate_alphabets(k)

                        to_automata.append([right + current, k + "_reverse", left + current])
                        backbone[k + "_reverse"] = reverse_alphabets(backbone[k])

        else:
            
            # to treat final endgroup
            to_automata.append([terminals[i - 1][1] + prev, end_group[i], right_desc + origin]) 
            backbone[end_group[i]] = generate_alphabets(end_group[i])
            backbone[end_group[i] + "_reverse"] = reverse_alphabets(backbone[end_group[i]])

        object_id_highest += 1

    for i in to_automata:
        fragments.append(i)
    
    print("fragments so far: ", fragments)

    if "start" not in to_automata[0][0]:

        # scale all indices by the object_id_highest

        o_id = 0
        prev = int(to_automata[0][0][4:])
        for i in range(len(to_automata)):

            scale0 = object_id_highest - prev + int(to_automata[i][0][4:])
            scale2 = object_id_highest - prev + int(to_automata[i][2][4:])

            o_id = max(o_id, int(to_automata[i][0][4:]) + scale0)
            o_id = max(o_id, int(to_automata[i][2][4:]) + scale2)

            scale0 = to_automata[i][0][0:4] + str(int(to_automata[i][0][4:]) + scale0)
            scale2 = to_automata[i][2][0:4] + str(int(to_automata[i][2][4:]) + scale2)

            if i == 0: # end group
                fragments.append([scale2, to_automata[i][1]+"_reverse", to_automata[i][0]])
            elif i == len(to_automata) - 1: # end group
                fragments.append([to_automata[i][2], to_automata[i][1]+"_reverse", scale0])
            else:
                fragments.append([scale0, to_automata[i][1], scale2])

        object_id_highest = o_id + 1
    
    print("after processing: ", fragments, object_id_highest)

    return fragments, backbone, object_id_highest

def reverse_paths(automata):
    a = []
    for i in automata:
        a.append(i[::-1])
    a = a[::-1]
    for i in range(len(a)):
        for j in range(len(a[i])):
            if a[i][j] == "start0":
                a[i][j] = "END"
            if a[i][j] == "end0":
                a[i][j] = "START"
    for i in range(len(a)):
        for j in range(len(a[i])):
            if a[i][j] == "START":
                a[i][j] = "start0"
            if a[i][j] == "END":
                a[i][j] = "end0"
    for i in range(len(a)):
        if "reverse" in a[i][1]:
            a[i][1] = a[i][1].replace("_reverse", "")
        elif a[i][1] != "[Cf]":
            a[i][1] = a[i][1] + "_reverse"
    return a

def generate_transitions(automata, backbone):
    polymer_graph = nx.DiGraph()

    for node in automata:
        compatible = get_compatible_descriptor(node[0])
        index_to_string = nx.get_node_attributes(polymer_graph, 'string')
        val_list = list(index_to_string.values())
        reverse = {value : key for (key, value) in index_to_string.items()}
        if compatible in val_list: 
            # if the compatible descriptor exists in the graph, only add SMILES
            new_node2 = polymer_graph.number_of_nodes()
            polymer_graph.add_node(new_node2, string = node[1], descriptor = False)
            polymer_graph.add_edge(reverse[compatible], new_node2)
        else:
            # if compatible does not exist, add it and the SMILES
            new_node1 = polymer_graph.number_of_nodes()
            polymer_graph.add_node(new_node1, string = compatible, descriptor = True)
            new_node2 = polymer_graph.number_of_nodes()
            polymer_graph.add_node(new_node2, string = node[1], descriptor = False)
            polymer_graph.add_edge(new_node1, new_node2)
        index_to_string = nx.get_node_attributes(polymer_graph, 'string')
        val_list = list(index_to_string.values())
        reverse = {value : key for (key, value) in index_to_string.items()}
        if node[2] in val_list:
            polymer_graph.add_edge(new_node2, reverse[node[2]])
        else:
            new_node3 = polymer_graph.number_of_nodes()
            polymer_graph.add_node(new_node3, string = node[2], descriptor = True)
            polymer_graph.add_edge(new_node2, new_node3)

    labels = nx.get_node_attributes(polymer_graph, 'string')
    initialState = ""
    finalState = ""
    exclude = -1

    # adjust the initial and final states if stochastic objects are fragments 
    for l in labels.items():
        if l[1] == "start0":
            initialState = l[0]
            for a in automata:
                if a[0] == "start0":
                    if a[1] == "[Cf]":
                        find = a[2]
                        for i in labels.items():
                            if i[1] == find:
                                initialState = i[0]
                                break
        if l[1] == "end0":
            finalState = l[0]
            for a in automata:
                if a[2] == "end0":
                    if a[1] == "[Cf]":
                        find = get_compatible_descriptor(a[0])
                        end = finalState
                        for i in labels.items(): 
                            neighbors = list(polymer_graph[i[0]])
                            if len(neighbors) > 0 and neighbors[0] == end: # eliminate [Cf] from graph
                                exclude = i[0]
                            if i[1] == find:
                                finalState = i[0]

    # visualize polymer graph
    # visual = nx.nx_pydot.write_dot(polymer_graph,"dot.txt")
    # text_file = open("dot.txt", "r")
    # data = text_file.read()
    # data = data.replace("string", "shape=circle, label")
    # graphs = pydot.graph_from_dot_data(data)
    # graphs[0].write_svg("output.svg")

    # compute transition function
    d = nx.get_node_attributes(polymer_graph, 'descriptor')
    n = nx.to_numpy_array(polymer_graph)
    n = np.argwhere(n == 1)
    adjacent_SMILES = []    
    for i in range(len(n)):
        if exclude in n[i]:
            continue
        one = n[i][0] # start index node
        pInitial = list(nx.all_simple_paths(polymer_graph, initialState, one))
        if len(pInitial) == 0 and one != initialState:
            continue
        if d[one]:
            two = n[i][1] #end index node
            for j in range(len(n)):
                if two == n[j][0]:
                    adjacent_SMILES.append([n[i], n[j]]) # in between SMILES

    statesList = []
    alphabet = []
    transition = []
    counter = nx.number_of_nodes(polymer_graph) # this is the first state for bb unit 1  
    for p in adjacent_SMILES:
        prev_state = p[0][0]
        smiles = p[0][1]
        next_state = p[1][1]
        bb = backbone[labels[smiles]]
        for s in range(len(bb)):
            bb[s] = Chem.MolToSmiles(Chem.MolFromSmiles(bb[s]))
            if bb[s] not in alphabet:
                alphabet.append(bb[s])
            if len(bb) == 1:
                step = [prev_state,next_state,bb[s]]
            else:
                counter_prev = counter
                counter += 1
                if s == 0:
                    step = [prev_state, counter, bb[s]]
                elif s == len(bb) - 1:
                    step = [counter_prev, next_state, bb[s]]
                else:
                    step = [counter_prev, counter, bb[s]]
            if str(step[0]) not in statesList:
                statesList.append(str(step[0]))
            if str(step[1]) not in statesList:
                statesList.append(str(step[1]))
            if [str(step[0]), str(step[1]), str(step[2])] not in transition:
                transition.append([str(step[0]), str(step[1]), str(step[2])])
    
    alphabet_dictionary = {}
    letter = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    counter = 0
    for a in alphabet:
        alphabet_dictionary[a] = letter[counter]
        counter += 1
    for t in range(len(transition)):
        transition[t][2] = alphabet_dictionary[transition[t][2]]

    initialState = str(initialState)
    finalState = str(finalState)
    alphabet = list(alphabet_dictionary.values())
    return initialState, finalState, statesList, alphabet, transition, alphabet_dictionary

if __name__ == "__main__":
    bigsmiles = "{[][>1]NCCN[>1],[<1]{[<1]O=C([>1])Nc1ccc(Cc2ccc(NC(=O)[>1])cc2)cc1,[<1]O{[<1]CC(C[>1])O[<1][>1]}[<1][<1]}[<1][]}"
    # bigsmiles = "OCCO{[>1][<1]CCO[>1],[<1]CC(C)O[>1][<1]}"
    # bigsmiles = "{[][$1]CC(C)[$1][]}"
    bigsmiles = process(bigsmiles)
    automata, backbone, o = generate_paths(bigsmiles, "start", "end", [], dict(), 1)
    print(automata)
    initialState, finalState, states, alphabet, transition_input, alphabet_dictionary = generate_transitions(automata, backbone)
    # print("initial, final: ",initialState, finalState)
    # print("states: ", states)
    # print("alphabet: ", alphabet)
    # print("transition_input: ", transition_input)
    # print("alphabet_dictionary: ", alphabet_dictionary)