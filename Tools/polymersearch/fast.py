import rdkit
from rdkit import Chem
import networkx as nx
import re

from .graphs import get_objects, get_repeats, desc_regex
from .search_tools import generate_NetworkX_graphs, identify_cycles, desc_regex

def generate_filters(bigsmiles):

    graphs = generate_NetworkX_graphs(bigsmiles, is_bigsmarts=False)
    group_cycles, individual_clusters, treat_as_cycle = identify_cycles(graphs)
    num_objects = len(group_cycles)

    contains_nested = False
    level = nx.get_node_attributes(graphs["atomistic"], "level")
    for l in level:
        if level[l] > 0:
            contains_nested = True
    
    num_desc = 0
    smarts = nx.get_node_attributes(graphs["multidigraph"], "symbol")
    for s in smarts:
        d = list(re.finditer(desc_regex, smarts[s]))
        if len(d) == 0:
            count = 0
            edges = graphs["multidigraph"].edges
            for e in edges:
                if s == e[0] or s == e[1]:
                    count += 1                
            if count > num_desc:
                num_desc = count

    rings = to_smiles(graphs) 

    return {
            "bigsmiles": bigsmiles, 
            "contains_nested": contains_nested, 
            "num_objects": num_objects, 
            "num_desc": num_desc, 
            "ring_smiles": rings,
            "atomistic": graphs["atomistic"],
            "topology": graphs["topology"],
            "top_undir": graphs["top_undir"],
            "multidigraph": graphs["multidigraph"]
            }

def cycle_to_ring_SMILES(graphs, cycle):
    # https://github.com/maxhodak/keras-molecules/pull/32/files     

    symbols = nx.get_node_attributes(graphs["atomistic"], "symbol")
    formal_charge = nx.get_node_attributes(graphs["atomistic"], "formal_charge")
    is_aromatic = nx.get_node_attributes(graphs["atomistic"], "is_aromatic")
    bonds = nx.get_edge_attributes(graphs["atomistic"], "bond_type_object")
    bonds_str = nx.get_edge_attributes(graphs["atomistic"], "bond_type")
    ids = nx.get_node_attributes(graphs["atomistic"], "ids")
    
    node_to_idx = dict()
    mol = Chem.RWMol()
    for node in graphs["atomistic"].nodes():
        if ids[node] in cycle and node not in graphs["descriptors"] and symbols[node] != "":
            a = Chem.Atom(symbols[node])
            a.SetFormalCharge(formal_charge[node])
            a.SetIsAromatic(is_aromatic[node]) 
            idx = mol.AddAtom(a)
            node_to_idx[node] = idx
    
    bonds_across_desc_obj = dict()
    already_added = set()
    for edge in graphs["atomistic"].edges():
        first, second = edge
        if ids[first] in cycle and ids[second] in cycle:
            if first in graphs["descriptors"] or second in graphs["descriptors"]:
                bonds_across_desc_obj[(first, second)] = bonds_str[tuple(sorted([first, second]))]
            else:
                ifirst = node_to_idx[first]
                isecond = node_to_idx[second]
                bond_type_object = bonds[first, second]
                if tuple(sorted([ifirst, isecond])) not in already_added:
                    mol.AddBond(ifirst, isecond, bond_type_object)
                    already_added.add(tuple(sorted([ifirst, isecond])))

    # already_added_descriptor = set()
    # for key1 in bonds_across_desc_obj:
    #     for key2 in bonds_across_desc_obj:
    #         if key1 == key2:
    #             continue
    #         descriptor = set(key1).intersection(key2)
    #         if len(descriptor) == 1: 
    #             try:
    #                 a = "1" in bonds_across_desc_obj[key1] and "2" in bonds_across_desc_obj[key2]
    #                 b = "2" in bonds_across_desc_obj[key1] and "1" in bonds_across_desc_obj[key2]
    #                 if not (a or b):
    #                     continue

    #                 a = set(key1).difference(descriptor).pop() 
    #                 b = set(key2).difference(descriptor).pop() 
    #                 descriptor = list(descriptor)[0] 
    #                 ifirst = node_to_idx[a] 
    #                 isecond = node_to_idx[b] 
    #                 if tuple(sorted([ifirst, isecond])) not in already_added and descriptor not in already_added_descriptor:
    #                     mol.AddBond(ifirst, isecond, rdkit.Chem.rdchem.BondType.SINGLE)
    #                     already_added.add(tuple(sorted([ifirst, isecond])))
    #                     already_added_descriptor.add(descriptor)
    #             except:
    #                 continue
    
    Chem.SanitizeMol(mol)
    mol = Chem.MolToSmiles(mol)
    return mol

def to_smiles(graphs):
    cycles, clusters, treat_as_cycle = identify_cycles(graphs)

    repeats = []
    egs = []

    # get all topological IDs
    ids = nx.get_node_attributes(graphs["atomistic"], "ids")
    ids = set(list(ids.values()))
    for id in ids:
        eg = True
        for object in clusters:
            if id in object:
                eg = False
                break
        if eg:
            smiles = cycle_to_ring_SMILES(graphs, [id])
            egs.append(smiles) 
            
    for stochastic_object in cycles: 
        smiles_in_object = []
        for cycle in stochastic_object:
            smiles = cycle_to_ring_SMILES(graphs, cycle)
            smiles_in_object.append(smiles)
        repeats.append(smiles_in_object)
    
    return repeats, egs

def search_small_molecules(graphs_q, graphs_t):
      
    # Do not canonicalize BigSMARTS
    atomistic = graphs_q["atomistic"]
    symbols = nx.get_node_attributes(atomistic, "symbol")
    for s in symbols: 
        if symbols[s] == "?*":
            symbols[s] = "Fm"
    nx.set_node_attributes(atomistic, symbols, "symbol")
    graphs_q[atomistic] = atomistic
    repeats_q, egs_q = to_smiles(graphs_q)
    repeats_t, egs_t = to_smiles(graphs_t)
    print(repeats_q)
    print(repeats_t)
    quit()

    # Linear + Concatenation
    # Branched
    # Localization
    # ?*

    def localization():
        atomistic = graphs_q["atomistic"]
        local_el = nx.get_node_attributes(atomistic, "local_el")
        for key in local_el:
            smarts_list = local_el[key]["ru_local_el"] 
            for smarts in smarts_list:
                q = Chem.MolFromSmarts(smarts)
                found = False
                for object_t in repeats_t:
                    for t in object_t:
                        t = Chem.MolFromSmiles(t)
                        matches = t.GetSubstructMatches(q)
                        if matches:
                            found = True
                            break
                if not found:
                    return False
        return True
            
    for object_q in repeats_q:
        for q in object_q:
            if "*[Fm]" == q: # wildcard object
                smarts = localization()
                if not smarts:
                    return False
            elif "[Fm]" in q: # ?* in RU
                print("test") 
            else: # repeat units
                q_obj = Chem.MolFromSmarts(q)
                found = False
                for object_t in repeats_t:
                    for t in object_t:
                        t_obj = Chem.MolFromSmiles(t)
                        matches = t_obj.GetSubstructMatches(q_obj)
                        if matches:
                            found = True
                            break
                if not found:
                    return False
            
    return True