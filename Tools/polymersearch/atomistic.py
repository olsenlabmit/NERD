import sys
import rdkit
import networkx as nx

from .adjacency import *

class Atom_Graph_Matcher:
    def __init__(self, G1, G2, all_descriptors_1, all_descriptors_2):
        self.G1 = G1
        self.G2 = G2
        self.G1_nodes = set(G1.nodes())
        self.G2_nodes = set(G2.nodes())
        self.G1_node_order = {n: i for i, n in enumerate(G1)}

        self.node_info_1 = (
                            nx.get_node_attributes(self.G1, "symbol"),
                            nx.get_node_attributes(self.G1, "formal_charge"),
                            nx.get_node_attributes(self.G1, "is_aromatic"),
                            nx.get_node_attributes(self.G1, "chirality"),
                            nx.get_node_attributes(self.G1, "num_hs"),
                            )
        self.node_info_2 = (
                            nx.get_node_attributes(self.G2, "symbol"),
                            nx.get_node_attributes(self.G2, "formal_charge"),
                            nx.get_node_attributes(self.G2, "is_aromatic"),
                            nx.get_node_attributes(self.G2, "chirality"),
                            nx.get_node_attributes(self.G2, "num_hs")
                            ) 

        self.bonds_1 = nx.get_edge_attributes(self.G1, "bond_type")
        self.bonds_2 = nx.get_edge_attributes(self.G2, "bond_type")
        self.ids_1 = nx.get_node_attributes(self.G1, "ids")
        self.ids_2 = nx.get_node_attributes(self.G2, "ids")
        self.descriptors_1 = all_descriptors_1
        self.descriptors_2 = all_descriptors_2
        self.q_macrocycle = []
        self.t_macrocycle = []

        # Set recursion limit.
        self.old_recursion_limit = sys.getrecursionlimit()
        expected_max_recursion_level = len(self.G1)
        if self.old_recursion_limit < 1.5 * expected_max_recursion_level:
            # Give some breathing room.
            sys.setrecursionlimit(int(1.5 * expected_max_recursion_level))

    def set_clusters(self, query_clusters, target_clusters, q_macrocycle, t_macrocycle):
        self.query_clusters = query_clusters
        self.target_clusters = target_clusters
        self.q_macrocycle = q_macrocycle
        self.t_macrocycle = t_macrocycle

    def init_dfs(self, cycle, cluster):
        self.partial_1 = {}
        self.partial_2 = {}
        self.candidates_1 = {}
        self.candidates_2 = {}
        self.state = GMState(self)
        if len(cycle) == 1:
            found = False
            for q_c in self.q_macrocycle:
                if list(cycle)[0] in q_c:
                    self.smarts_endgroup == False
                    cycle = q_c
                    found = True
            if not found:
                self.smarts_endgroup = True
        else:
            self.smarts_endgroup = False

        # cycle and cluster information for localization search
        self.cycle = cycle
        self.cluster = cluster
        self.original_start = -1
        self.end1 = False
        self.end2 = False

        # initialize information for traversing across descriptors
        self.neighbors_1 = dict()
        self.neighbors_2 = dict()
        self.traversal_2 = dict()
        for key in self.ids_2:
            self.traversal_2[self.ids_2[key]] = 0
        
        # get total number of query atoms in the cycle as a stopping point
        self.query_atoms = 0 
        for key in self.node_info_1[0]:
            if self.ids_1[key] in self.cycle and key not in self.descriptors_1:
                self.query_atoms += 1
        
    def generate_matches(self):
        yield from self.level()
        
    def level(self):
        if len(self.partial_1) == self.query_atoms:

            def target_cycle_has_descriptor(processed):

                has_wildcard = False
                for key in self.partial_1:
                    if type(self.partial_1[key]) == list:
                        has_wildcard = True
                        break
                
                if has_wildcard:
                    for key in self.partial_1:
                        match = self.partial_1[key]
                        if type(match) != list:
                            match = [match]
                        for key2 in match:
                            if key2 in self.descriptors_2:
                                return True
                    fragment_set = set()
                    for key in self.partial_1:
                        if type(self.partial_1[key]) != list:
                            fragment_set.add(get_iteration(self.partial_1[key]))
                        else:
                            for key in self.partial_1[key]:
                                fragment_set.add(get_iteration(key))
                        if len(fragment_set) >= 2:
                            return True
                    return False

                q_ends_con = directly_connected(self.start_atom, self.end_atom, self.G1)
                q_ends_sep_desc = descriptors_separating(self.start_atom, self.end_atom, self.bonds_1, self.descriptors_1)
                
                start_atom_mapped = processed[self.start_atom]
                end_atom_mapped = processed[self.end_atom]
                t_ends_con = directly_connected(start_atom_mapped, end_atom_mapped, self.G2)
                t_ends_sep_desc = descriptors_separating(start_atom_mapped, end_atom_mapped, self.bonds_2, self.descriptors_2)

                if len(q_ends_sep_desc) >= 2 or q_ends_con and len(q_ends_sep_desc) >= 1:
                    if len(t_ends_sep_desc) >= 2 or t_ends_con and len(t_ends_sep_desc) >= 1:
                        return True
            
                same = True
                fragment_set = set()
                for key in self.partial_1:
                    if type(self.partial_1[key]) != list:
                        fragment_set.add(get_iteration(self.partial_1[key]))
                        if len(fragment_set) >= 2:
                            same = False
                            break

                if same: 
                    if len(t_ends_sep_desc) >= 1:
                        return True
                else: 
                    if t_ends_con or len(t_ends_sep_desc) >= 1:
                        return True

                return False
            
            def bond_adjacency(processed):
                if self.original_start != -1 or self.node_info_1[0][self.end_atom] == "?*":
                    return True
                if "=" in self.node_info_1[0][self.start_descriptor]:
                    a = processed[self.start_atom]
                    b = processed[self.end_atom]
                    desc = descriptors_separating(a, b, self.bonds_2, self.descriptors_2)
                    if desc:
                        if "=" in self.node_info_2[0][desc[0]]:
                            return True
                        return False
                    if "DOUBLE" in self.bonds_2[feed_to_bonds_n(a, b)]:
                        return True
                    return False
                else:
                    a = processed[self.start_atom]
                    b = processed[self.end_atom]
                    desc = descriptors_separating(a, b, self.bonds_2, self.descriptors_2)
                    if desc:
                        if "=" in self.node_info_2[0][desc[0]]:
                            return False
                        return True
                    if "DOUBLE" in self.bonds_2[feed_to_bonds_n(a, b)]:
                        return False
                    return True

            processed = dict()
            for key in self.partial_1:
                if type(self.partial_1[key]) == str:
                    processed[key] = get_node_index(self.partial_1[key])
                else:
                    processed[key] = [] 
                    for mapping in self.partial_1[key]:
                        if type(mapping) == int:
                            processed[key].append(mapping)
                        else:
                            processed[key].append(get_node_index(mapping))
            if self.smarts_endgroup:
                yield processed 
            else:
                if enforce_adjacency(processed, self, self.cycle) and bond_adjacency(processed) and target_cycle_has_descriptor(processed):
                    yield processed 

        else:
            for G1_node, G2_node, path in self.suggest_nodes():
                if self.do_nodes_match(G1_node, G2_node, path):
                    # print("Pass Semantic Checks: ", G1_node, G2_node, path)
                    newstate = self.state.__class__(self, G1_node, G2_node, path)
                    # print("Partial solution: ", self.partial_1)
                    yield from self.level()
                    newstate.backtrack()
                        
    def suggest_nodes(self):
        # generate query and target candidate atoms that are not in the partial solution
        candidates_query = [node for node in self.candidates_1 if node not in self.partial_1]
        candidates_target = [node for node in self.candidates_2 if node not in self.partial_2]

        # delete first atom if it is ?*
        if self.original_start in candidates_query:
            candidates_query.remove(self.original_start)

        if len(self.partial_1) == self.query_atoms - 1 and self.original_start != -1:
            yield from self.suggest_paths()

        # if both are populated, then pick a query node and suggest target nodes
        if candidates_query and candidates_target:
            min_key = self.G1_node_order.__getitem__
            node_1 = min(candidates_query, key=min_key)
            if self.node_info_1[0][node_1] == "?*":
                yield from self.suggest_paths(node_1)
            else:
                for node_2 in candidates_target:
                    yield node_1, node_2, None
        
        else:
            if self.smarts_endgroup:
                for key in self.node_info_1[0]:
                    if self.ids_1[key] in self.cycle:
                        self.start_atom = key
                        self.start_descriptor = -1
                        self.end_atom = key
                        break
            else:
                found = False
                for key in self.node_info_1[0]:
                    if self.ids_1[key] in self.cycle and key in self.descriptors_1:
                        neighbors = self.G1[key]                        
                        for start in neighbors:
                            for end in neighbors:
                                desc = descriptors_separating(start, end, self.bonds_1, self.descriptors_1)
                                if key in desc and start != end and self.ids_1[start] in self.cycle and self.ids_1[end] in self.cycle and not found:
                                    self.start_atom = start
                                    self.start_descriptor = key
                                    self.end_atom = end
                                    found = True

                if self.node_info_1[0][self.start_atom] == "?*":
                    neighbors = self.G1[self.start_atom]     
                    for new_start in neighbors:
                        if new_start not in self.descriptors_1:
                            self.original_start = self.start_atom
                            self.start_atom = new_start
                            break
                # print("Starting positions: ", self.start_atom, self.start_descriptor, self.end_atom)
                    
            # once the query atoms are chosen, choose the target atom
            for target_start in self.node_info_2[0]:
                if target_start not in self.partial_2 and self.ids_2[target_start] in self.cluster and target_start not in self.descriptors_2:
                    if "_" not in str(target_start):
                        target_start = str(target_start) + "_0"
                    yield self.start_atom, target_start, None

    def suggest_paths(self, wildcard = -1):

        def path_valid(path):
            for i in range(len(path)):
                if i > 0 and i < len(path) - 1 and path[i] in self.descriptors_2:
                    pair_1 = self.bonds_2[feed_to_bonds_n(path[i - 1], path[i])]
                    pair_2 = self.bonds_2[feed_to_bonds_n(path[i], path[i + 1])]
                    if pair_1 != get_comp_descriptor(pair_2):
                        return False
            return True
        
        def add_iteration(path):
            desc_crossed = False
            for i in range(len(path)):
                if path[i] in self.descriptors_2:
                    desc_crossed = True
                    continue
                if path[i] not in self.descriptors_2 and not desc_crossed:
                    path[i] = str(path[i]) + "_" + get_iteration(self.partial_1[back_q])
                elif path[i] not in self.descriptors_2 and desc_crossed:
                    path[i] = str(path[i]) + "_" + str(self.traversal_2[self.ids_2[path[i]]] + 1)
            return path

        def exist_in_partial(path, forward_t):
            path = path + [forward_t]
            for i in range(len(path)):
                if path[i] not in self.descriptors_2:
                    for key in self.partial_1:
                        already = self.partial_1[key]
                        if type(already) == str:
                            already = [already]
                        if path[i] in already:
                            return True
            return False

        if wildcard == -1:
            if self.node_info_1[0][self.end_atom] == "?*":
                forward_q = self.start_atom
                forward_t = self.partial_1[self.start_atom]
                self.end2 = True
                yield forward_q, forward_t, (self.original_start, self.partial_1[self.end_atom])
            else:
                forward_q = self.start_atom
                forward_t = self.partial_1[self.start_atom]
                forward_t = get_node_index(forward_t)
                back_q = self.end_atom
                back_t = self.partial_1[self.end_atom]
                back_t = get_node_index(back_t)
                for path in nx.all_simple_paths(self.G2, source = back_t, target = forward_t):
                    if path_valid(path):
                        path = add_iteration(path)
                        forward_t = path[-1]
                        path = path[1:-1]
                        if path == []:
                            continue
                        if not exist_in_partial(path, forward_t):
                            self.end1 = True
                            yield forward_q, forward_t, (self.original_start, path)
            
        elif wildcard == self.end_atom:
            neighbors = self.G1[wildcard]
            partial_solution = list(self.partial_1.keys()) 

            forward_q = self.start_atom
            forward_t = self.partial_1[self.start_atom]
            forward_t = get_node_index(forward_t)
            for node in neighbors:
                if node in partial_solution:
                    back_q = node 
            back_t = self.partial_1[back_q]
            back_t = get_node_index(back_t)

            if forward_t == back_t:
                neighbors = self.G2[back_t]
                for n in neighbors:
                    if neighbors not in self.descriptors_2:
                        back_t = n
                        break
            for path in nx.all_simple_paths(self.G2, source = back_t, target = forward_t):
                if path_valid(path):
                    path = add_iteration(path)
                    forward_t = path[-1]
                    path = path[1:-1]
                    if path == []:
                        continue
                    if not exist_in_partial(path, forward_t):
                        self.end1 = True
                        yield forward_q, forward_t, (self.end_atom, path)
        
        else:
            neighbors = self.G1[wildcard]
            partial_solution = list(self.partial_1.keys()) 
            for node in neighbors:
                if node in partial_solution:
                    back_q = node 
                else:
                    forward_q = node
            back_t = self.partial_1[back_q]
            back_t = get_node_index(back_t)

            forward_t_candidates = nx.bfs_edges(self.G2, source = back_t)
            forward_t_candidates = set([i for sub in forward_t_candidates for i in sub])
            for forward_t_candidate in forward_t_candidates:
                a = forward_q in self.descriptors_1 and forward_t_candidate in self.descriptors_2
                b = forward_q not in self.descriptors_1 and forward_t_candidate not in self.descriptors_2

                if a or b: 
                    for path in nx.all_simple_paths(self.G2, source = back_t, target = forward_t_candidate):
                        if path_valid(path):
                            path = add_iteration(path)
                            forward_t = path[-1]
                            path = path[1:-1]
                            if not exist_in_partial(path, forward_t):
                                yield forward_q, forward_t, (wildcard, path)

    def do_nodes_match(self, G1_candidate, G2_candidate, path):
        G2_candidate_node_index = get_node_index(G2_candidate)
        G2_candidate_node_iteration = get_iteration(G2_candidate)
        
        # check whether atom symbol, formal charge, aromaticity, and chirality match
        for i in range(3):
            if self.node_info_1[0][G1_candidate] != "*": 
                if self.node_info_1[i][G1_candidate] != self.node_info_2[i][G2_candidate_node_index]: 
                    return False    
        if self.node_info_1[3][G1_candidate] != rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            if self.node_info_1[3][G1_candidate] != self.node_info_2[3][G2_candidate_node_index]: 
                return False
        # check whether hydrogens match
        if self.node_info_1[4][G1_candidate] > self.node_info_2[4][G2_candidate_node_index]:
            return False
        
        if path != None:
            return True

        # Check if neighboring atoms and bonds are the same
        # get previous neighbor and descriptor 
        G1_neighbor_before_descriptor = None
        desc_id_1 = None
        for key in self.neighbors_1:
            if G1_candidate in self.neighbors_1[key]:
                # are there any atoms in the same fragment that are part of the same solution
                add = True
                for k in self.neighbors_1[key]:
                    if self.ids_1[k] == self.ids_1[G1_candidate] and k != G1_candidate and k in list(self.partial_1.keys()):
                        add = False
                if add:
                    G1_neighbor_before_descriptor = key[0]
                    desc_id_1 = key[1]

        # get previous neighbor and descriptor 
        G2_neighbor_before_descriptor = None
        desc_id_2 = None
        for key in self.neighbors_2:
            if G2_candidate in self.neighbors_2[key]:
                # are there any atoms in the same fragment that are part of the same solution
                add = True
                for k in self.neighbors_2[key]:
                    if self.ids_2[get_node_index(k)] == self.ids_2[G2_candidate_node_index] and get_node_index(k) != G2_candidate_node_index and k in self.partial_1:
                        add = False
                if add:
                    G2_neighbor_before_descriptor = key[0]
                    desc_id_2 = key[1]

        for neighbor in self.G1[G1_candidate]:
            # when getting neighbors in query graph, if neighbor is descriptor, get atom before it
            if neighbor == desc_id_1:
                pair_1 = self.bonds_1[feed_to_bonds_n(G1_candidate, desc_id_1)][2:]
                neighbor = G1_neighbor_before_descriptor

            elif neighbor in self.partial_1:
                pair_1 = self.bonds_1[feed_to_bonds_n(G1_candidate, neighbor)]
            
            # only check bond connections in the core
            if neighbor in self.partial_1:

                mapped = self.partial_1[neighbor]
                mapped_node_index = get_node_index(mapped)
                mapped_node_iteration = get_iteration(mapped)

                if mapped == G2_neighbor_before_descriptor:
                    pair_2 = self.bonds_2[feed_to_bonds_n(mapped_node_index, desc_id_2)][2:]
                else:
                    if G2_candidate_node_iteration == mapped_node_iteration:
                        try:
                            pair_2 = self.bonds_2[feed_to_bonds_n(mapped_node_index, G2_candidate_node_index)]
                        except:
                            return False
                    else:
                        return False

                if pair_1 != pair_2:
                    return False
                
        # Following a descriptor node, only one fragment is allowed to be traversed. 
        for key in self.neighbors_2:
            if G2_candidate in self.neighbors_2[key]:
                for node in self.neighbors_2[key]:
                    if node != G2_candidate and node in self.partial_2:
                        return False

        return True
        
class GMState:  
    def __init__(self, GM, G1_candidate = None, G2_candidate = None, path = None):
        self.GM = GM

        # Initialize the last stored node pair.
        self.G1_node = None
        self.G2_node = None
        self.depth = len(GM.partial_1) 

        if G1_candidate is None or G2_candidate is None:
            GM.partial_1 = {}
            GM.partial_2 = {}
            GM.candidates_1 = {}
            GM.candidates_2 = {}

        if G1_candidate is not None and G2_candidate is not None:      
            # add candidates to partial solution
            if G1_candidate not in GM.partial_1 and G1_candidate not in GM.descriptors_1:
                GM.partial_1[G1_candidate] = G2_candidate
            if G2_candidate not in GM.partial_2 and G2_candidate not in GM.descriptors_2:
                GM.partial_2[G2_candidate] = G1_candidate
            if path != None and path[0] not in GM.partial_1:
                GM.partial_1[path[0]] = path[1]

            # store new candidates for backtracking purposes
            self.G1_node = G1_candidate 
            self.G2_node = G2_candidate 
            self.path = path

            # augment candidates to select for DFS
            self.depth = len(GM.partial_1)
            if G1_candidate not in GM.candidates_1 and G1_candidate not in GM.descriptors_1:
                GM.candidates_1[G1_candidate] = self.depth
            if G2_candidate not in GM.candidates_2 and G2_candidate not in GM.descriptors_2:
                GM.candidates_2[G2_candidate] = self.depth
            
            # update neighbors
            for node in GM.partial_1:
                if type(GM.partial_1[node]) == list:
                    continue
                for neighbor in GM.G1[node]:
                    if neighbor in GM.descriptors_1: 
                        if neighbor != GM.start_descriptor:
                            # if the node and its neighbor are already in neighbors_1, do not add to candidates
                            cont = False
                            for key in GM.neighbors_1:
                                if (node, neighbor) == key or node in GM.neighbors_1[key] and neighbor == key[1]:
                                    cont = True
                                    break
                            if cont:
                                continue
                            
                            # otherwise add to neighbors_1
                            GM.neighbors_1[(node, neighbor)] = []
                            pair = feed_to_bonds_n(node, neighbor)
                            for bond_pair in GM.bonds_1:
                                if neighbor in bond_pair and GM.bonds_1[bond_pair] == get_comp_descriptor(GM.bonds_1[pair]):
                                    new_node = [i for i in bond_pair if i != neighbor][0] 
                                    if GM.ids_1[new_node] in GM.cycle:
                                        if new_node not in GM.candidates_1:
                                            GM.candidates_1[new_node] = self.depth
                                        
                                        # GM neighbors maps node and neighbor to list of nodes across the descriptor
                                        GM.neighbors_1[(node, neighbor)].append(new_node)
                    else:
                        # add to candidates if not already there and not a descriptor
                        if neighbor not in GM.candidates_1 and neighbor not in GM.partial_1 and neighbor not in GM.descriptors_1:
                            GM.candidates_1[neighbor] = self.depth

            for original in GM.partial_2:
                node_index = get_node_index(original)
                iteration = get_iteration(original)
                for neighbor in GM.G2[node_index]:
                    if neighbor in GM.descriptors_2:

                        # if the node and its neighbor are already in neighbors_2, do not add to candidates
                        cont = False
                        for key in GM.neighbors_2:
                            if (original, neighbor) == key or original in GM.neighbors_2[key] and neighbor == key[1]:
                                cont = True
                                break
                        if cont:
                            continue

                        # otherwise add to neighbors_2
                        GM.neighbors_2[(original, neighbor)] = []
                        pair = feed_to_bonds_n(node_index, neighbor)
                        for bond_pair in GM.bonds_2:
                            if neighbor in bond_pair and GM.bonds_2[bond_pair] == get_comp_descriptor(GM.bonds_2[pair]):
                                new_node = [i for i in bond_pair if i != neighbor][0]
                                id = GM.ids_2[new_node]
                                if id in GM.cluster:
                                    GM.traversal_2[id] += 1
                                    label = str(new_node) + "_" + str(GM.traversal_2[id])
                                    GM.candidates_2[label] = self.depth
                                    GM.neighbors_2[(original, neighbor)].append(label)
                    else:
                        # add to candidates if not already there and not a descriptor
                        label = str(neighbor) + "_" + iteration
                        if label not in GM.candidates_2 and label not in GM.partial_2 and neighbor not in GM.descriptors_2:
                            GM.candidates_2[label] = self.depth

    def backtrack(self):
        if self.G1_node is not None and self.G2_node is not None and not self.GM.end1 and not self.GM.end2:
            del self.GM.partial_1[self.G1_node]
            del self.GM.partial_2[self.G2_node]

        if self.path is not None:
            del self.GM.partial_1[self.path[0]]
            if self.GM.end2:
                self.GM.end2 = False
                self.GM.end1 = True
            elif self.GM.end1:
                self.GM.end1 = False

        for node in self.GM.candidates_1.copy():
            if self.GM.candidates_1[node] == self.depth:
                del self.GM.candidates_1[node]
                for key in self.GM.neighbors_1.copy():
                    if node in self.GM.neighbors_1[key]:
                        del self.GM.neighbors_1[key]

        for node in self.GM.candidates_2.copy():
            if self.GM.candidates_2[node] == self.depth:
                del self.GM.candidates_2[node]
                for key in self.GM.neighbors_2.copy():
                    if node in self.GM.neighbors_2[key]:
                        for val in self.GM.neighbors_2[key]:
                            node_index = get_node_index(val)
                            id = self.GM.ids_2[node_index] 
                            self.GM.traversal_2[id] -= 1
                        del self.GM.neighbors_2[key]