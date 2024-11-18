import rdkit
from rdkit import Chem
import networkx as nx
import re
import copy
import pandas as pd
import time

from .topology import Topology_Graph_Matcher
from .graphs import build_atomistic, build_topology, get_objects, get_repeats, desc_regex, ladder_desc_regex, atomistic_ladder
from .adjacency import feed_to_bonds_n

def generate_NetworkX_graphs(input_string, is_bigsmarts = False):
    ladder = list(re.finditer(ladder_desc_regex,input_string))
    if ladder:
        return generate_NetworkX_ladder_graphs(input_string, is_bigsmarts)
    atomistic_directed, no_other_objects = build_atomistic(input_string, is_bigsmarts)
    topology, topology_undir, multidigraph, descriptors, ids = build_topology(atomistic_directed)
    graphs = {"string": input_string, 
            "atomistic": atomistic_directed,
            "topology": topology, 
            "top_undir": topology_undir, 
            "ids": ids,
            "descriptors": descriptors, 
            "multidigraph": multidigraph,
            "no_other_objects":no_other_objects
            }
    return graphs

def generate_NetworkX_ladder_graphs(input_string, is_bigsmarts = False):
    atomistic = atomistic_ladder(input_string)
    atomistic_top = copy.deepcopy(atomistic)
    bond_type = nx.get_edge_attributes(atomistic_top, "bond_type")
    for b in bond_type:
        if bond_type[b] == "1a":
            bond_type[b] = "1"
        if bond_type[b] == "2a":
            bond_type[b] = "2"
        if bond_type[b] == "1b":
            bond_type[b] = "1"
        if bond_type[b] == "2b":
            bond_type[b] = "2"
    nx.set_edge_attributes(atomistic_top, bond_type, "bond_type")
    topology, topology_undir, multidigraph, descriptors, ids = build_topology(atomistic_top)
    ids = nx.get_node_attributes(atomistic_top, "ids")
    nx.set_node_attributes(atomistic, ids, "ids")
    graphs = { 
                "atomistic": atomistic,
                "topology": topology, 
                "top_undir": topology_undir, 
                "ids": ids, 
                "descriptors": descriptors, 
                "multidigraph": multidigraph
            }
    return graphs

def identify_cycles(graphs):
    topology = graphs["topology"]

    ids = nx.get_node_attributes(topology, "ids")
    individual_cycles = []
    for i in nx.simple_cycles(topology):
        individual_cycles.append([ids[j] for j in i])
    
    explicit = nx.get_node_attributes(graphs["atomistic"], "explicit_atom_ids")
    ids = nx.get_node_attributes(graphs["atomistic"], "ids")
    treat_as_cycle = []
    for key in explicit:
        if explicit[key]:
            remove = -1
            for i in individual_cycles:
                if ids[key] in i:
                    remove = i
                    break
            if remove != -1:
                individual_cycles.remove(remove)
                treat_as_cycle.append(remove)
    
    # https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements
    def to_graph(l):
        G = nx.Graph()
        for part in l:
            G.add_nodes_from(part)
            G.add_edges_from(to_edges(part))
        return G

    def to_edges(l):
        it = iter(l)
        last = next(it)
        for current in it:
            yield last, current
            last = current   
    
    individual_clusters = []
    G = to_graph(individual_cycles)
    for cycles in nx.connected_components(G):
        individual_clusters.append(cycles)
        
    level = nx.get_node_attributes(graphs["top_undir"], "level")
    ids = nx.get_node_attributes(graphs["top_undir"], "ids")
    inv_map = {v: k for k, v in ids.items()}
    bonds = nx.get_edge_attributes(graphs["top_undir"],"bond_type")
    nested = []
    for i in range(len(individual_clusters)):
        count_0 = 0
        count_1 = 0
        for j in individual_clusters[i]:
            id = inv_map[j]
            if id in level and level[id] == 1:
                count_1 = 1
            if id in level and level[id] == 0:
                count_0 = 1
        if count_0 == 0 and count_1 == 1 and individual_clusters[i] not in nested:
            nested.append(individual_clusters[i]) 
    for i in range(len(nested)):
        if nested[i] not in individual_clusters:
            continue
        for j in range(len(individual_clusters)):
            adjacent = False
            for k in nested[i]:
                for l in individual_clusters[j]:
                    if feed_to_bonds_n(inv_map[k], inv_map[l]) in bonds:
                        adjacent = True
            if adjacent:
                individual_clusters[j].update(nested[i])
                individual_clusters.remove(nested[i])
                break

    group_cycles = []
    for i in range(len(individual_clusters)):
        group_cycles.append([])
        for j in range(len(individual_cycles)):
            if individual_cycles[j][0] in individual_clusters[i]:
                group_cycles[-1].append(individual_cycles[j]) 
    return group_cycles, individual_clusters, treat_as_cycle

def contains_substructure(bigsmarts_graphs, bigsmiles_graphs):
    # the first step is identifying cycles or repeat units in the directed topology graphs for the query and target
    q_cycles, q_clusters, q_macrocycle = identify_cycles(bigsmarts_graphs)   
    # q_cycles is a list of list of lists (# objects x # cyles in object x # topology ID elements in that cycle)
    # q_clusters is a list of lists (# objects x # topology ID elements in that object)
    # nx.simple_cycles to detect cycles and nx.connected_components to identify adjacent cycles (stochastic objects)
    # q__macrocycle contains extra information for macrocycle polymer chemistries, in which the end group 
    # forms a cycle but is not treated as a repeat unit!
    t_cycles, t_clusters, t_macrocycle = identify_cycles(bigsmiles_graphs) 

    # the variable search is instantiated with all parameters needed to run graph traversal on the atomistic graph
    search = Topology_Graph_Matcher(bigsmarts_graphs, bigsmiles_graphs, 
                                            q_cycles, q_clusters, 
                                            t_cycles, t_clusters, 
                                            q_macrocycle, t_macrocycle)
    
    # the graph traversal is run here, which returns True or False
    return search.search_repeats_endgroups() 
    
def logical_repeat_unit_search(bigsmarts):
    # Assumptions: 
    # only 1 object in the query has logical operators
    # only one type of logical operation per query: "!", "or", "xor" along with "and"

    # determine all objects in the string
    objects = get_objects(bigsmarts)

    # iterate through each object
    for object in objects[0]:

        # get all repeat units in the object
        repeats = get_repeats(object)

        # map the logical operator to the repeat units
        logical = {} 
        for repeat in repeats[0]:

            # group logical operator with repeat units or SMARTS
            if repeat.find("[or") == 0 and repeat[0:5] not in logical:
                logical[repeat[0:5]] = [repeat[5:]]
            elif repeat.find("[or") == 0 and repeat[0:5] in logical:
                logical[repeat[0:5]].append(repeat[5:])
            
            elif repeat.find("[xor") == 0 and repeat[0:6] not in logical:
                logical[repeat[0:6]] = [repeat[6:]]
            elif repeat.find("[xor") == 0 and repeat[0:6] in logical:
                logical[repeat[0:6]].append(repeat[6:])

            elif repeat.find("!") == 0 and repeat != "!*" and "!" not in logical:
                logical["!"] = [repeat[1:]]
            elif repeat.find("!") == 0 and repeat != "!*" and repeat[1:] in logical:
                logical["!"].append(repeat[1:])

            elif "and" not in logical:
                logical["and"] = [repeat]
            else:
                logical["and"].append(repeat)
        
        # is this the object with the logical operators?
        logic = list(logical.keys())
        logic = [i for i in logic if i != "and"]

        # if not, continue
        if len(logic) == 0:
            continue
        
        # print("logical: ", logical)
        # list of object strings that convert logical strings into valid BigSMARTS
        objects_replaced = []
        logic_return = "and"
        for logic in logical:

            if "or" in logic or "xor" in logic:
                logic_return = logic
                for repeat in logical[logic]:
                    # delete logical operator and repeat unit
                    if logic + repeat + "," in object:
                        # this is for every other repeat unit in the stochastic object
                        replaced = object.replace(logic + repeat + ",", "")
                    else:
                        # this is for the last repeat unit in the stochastic object only
                        replaced = object.replace("," + logic + repeat, "")
                    replaced = replaced.replace(logic, "")
                    objects_replaced.append(replaced)
            
            elif "!" in logic:
                logic_return = logic
                for repeat in logical[logic]:
                    # delete logical operator
                    replaced = object.replace("!", "")
                    objects_replaced.append(replaced)
                # delete logical operator and repeat unit
                if logic + repeat + "," in object:
                    # this is for every other repeat unit in the stochastic object
                    replaced = object.replace(logic + repeat + ",", "")
                else:
                    # this is for the last repeat unit in the stochastic object only
                    replaced = object.replace("," + logic + repeat, "")
                replaced = replaced.replace(logic, "")
                objects_replaced.append(replaced)

        bigsmarts_replaced = []
        for o in objects_replaced:
            bigsmarts_replaced.append(bigsmarts.replace(object, o))
    
        return bigsmarts_replaced, logic_return
    
    return [bigsmarts], "and"
                
def graph_traversal_match(bigsmarts, bigsmiles):
    # this inputs a BigSMARTS and generates a set of equivalent strings that need to be searched
    # if the user types in the logical "or" or {[][or1][<]CCO[>],[<]CC(C)O[>][]}, 
    # then this function will generate two equivalent strings, one for PEG: {[][<]CCO[>][]} and one for PPO {[][<]CC(C)O[>][]} Either must be present!
    # if the user enters {[]![<]CCO[>][]}, then the function returns {[][<]CCO[>][]}, which is searched
    bigsmarts_list, logic = logical_repeat_unit_search(bigsmarts) 
    matches = []
    for bigsmarts in bigsmarts_list:
        # generate NetworkX molecular graphs, returns dictionary that includes atomistic and topology NetworkX graphs
        # this function should be called upon ingestion into CRIPT
        # if is_bigsmarts is True, then the graph generated has no information on hydrogens because molecular queries do not include hydrogens unless the user specified them
        # for example, the query CCO searches carbon-carbon-oxygen with single bond connections and no hydrogens
        # the target CCO searches the molecule ethanol or [CH3][CH2][OH]
        # atomistic undirected graph: each node has element symbol ("C" for carbon), each edge has bond order (single, double, triple)
        # topology directed graph: each SMARTS, SMILES, or bonding descriptor is compressed into a single node with a topological ID
        bigsmarts_graphs = generate_NetworkX_graphs(input_string = bigsmarts, is_bigsmarts = True) 
        bigsmiles_graphs = generate_NetworkX_graphs(input_string = bigsmiles, is_bigsmarts = False) 

        # runs graph traversal and returns True or False 
        m = contains_substructure(bigsmarts_graphs = bigsmarts_graphs, bigsmiles_graphs = bigsmiles_graphs)
        
        # the following statements are conditions in case the user specified logical queries
        if "or" in logic and "xor" not in logic and m:
            return True
        elif "and" in logic:
            return m
        matches.append(m)

    if "!" in logic:
        if matches[0] == False and matches[1] == False:
            return True
        return False 
    if "xor" in logic:
        if matches[0] == False and matches[1] == True or matches[0] == True and matches[1] == False:
            return True
        return False
    if "or" in logic:
        return matches[1]
    return False

def bigsmarts_unit_testing(bigsmarts_testing, sheet_name):
    data = pd.read_excel("polymersearch/paper/Validation.xlsx", sheet_name = sheet_name)
    data = data.fillna(0)
    if bigsmarts_testing == [0, -1]:
        queries = list(data.columns)[4:][0:]
    elif type(bigsmarts_testing) == list:
        queries = list(data.columns)[4:][bigsmarts_testing[0]:bigsmarts_testing[1]]
    else:
        queries = [bigsmarts_testing]
    targets = list(data["Targets (BigSMILES)"][3:])
    incorrect = False
    with open("Validation_" + sheet_name + ".txt", 'w') as f:
        tic = time.perf_counter()
        f.write("Start Time\n")
        # playsound("polymersearch/paper/Timer_End.mp3")
        for query in queries:
            f.write("Query: " + query + "\n")
            actual = list(data[query][3:])
            for i in range(len(targets)):
                if i % 50 == 0:
                    f.write("# of targets checked: " + str(i) + "/" + str(len(targets)) + "\n")
                try:
                    predicted = graph_traversal_match(query, targets[i])
                    if predicted != actual[i]:
                        f.write("INCORRECT: " + str(i + 5) + targets[i] + "\n")
                        incorrect = True
                except:
                    f.write("PROGRAM ERROR: " + str(i + 5) + targets[i] + "\n")
        toc = time.perf_counter()
        f.write("Stop Time\n")
        f.write(f"Total time is: {(toc - tic):0.4f} seconds\n")
        f.write(f"Average time is: {(toc - tic)/len(queries):0.4f} seconds\n")
        f.write("Is anything incorrect?", incorrect, "\n")
        # playsound("polymersearch/paper/Timer_End.mp3")