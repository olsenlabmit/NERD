import itertools
import copy

from .atomistic import Atom_Graph_Matcher
from .search_tools import *
from .adjacency import *

class Topology_Graph_Matcher:
    
    def __init__(self, q_graphs, t_graphs, q_cycles, q_clusters, t_cycles, t_clusters, q_macrocycle, t_macrocycle):        
        # here, the parameters for graph traversal are stored
        # input BigSMILES and BigSMARTS strings
        self.bigsmarts = q_graphs["string"]
        self.bigsmiles = t_graphs["string"]

        # query and target atomistic and toplogy graphs
        self.q_atom_graph = q_graphs["atomistic"]
        self.t_atom_graph = t_graphs["atomistic"]
        self.q_top_undir_graph = q_graphs["top_undir"]
        self.t_top_undir_graph = t_graphs["top_undir"]

        # a list of nodes in the atomistic graphs (atomg) that are bonding descriptors
        self.q_atomg_nodeind_desc = q_graphs["descriptors"]
        self.t_atomg_nodeind_desc = t_graphs["descriptors"]

        # BigSMARTS grammatically allows users to enter !{[][]}, which means no other objects in the target, this is a Boolean
        self.no_other_objects = q_graphs["no_other_objects"]

        # macrocycle information
        self.q_macrocycle = q_macrocycle
        self.t_macrocycle = t_macrocycle

        # define a variable for the graph traversal, where a depth-first search is run
        self.atomistic_dfs = Atom_Graph_Matcher(q_graphs["atomistic"], 
                                                        t_graphs["atomistic"], 
                                                        q_graphs["descriptors"], 
                                                        t_graphs["descriptors"]) 

        # maps node id to id that links atomistic to topology
        # 1 is query, 2 is target
        # these are dictionaries that map the atomistic ID to its topological ID
        # this is important because repeat unit cycles are searched first
        self.q_atom_ids = self.atomistic_dfs.ids_1 
        self.t_atom_ids = self.atomistic_dfs.ids_2
        self.all_ids_2 = list(self.t_atom_ids.values())
        
        # these are from the function identify_cycles
        self.query_cycles = q_cycles
        self.query_clusters = q_clusters
        self.target_cycles = t_cycles
        self.target_clusters = t_clusters

        # define the query and target end groups
        query_clusters = set()
        for i in self.query_clusters:
            query_clusters.update(i) 
        self.query_endgroups = set()
        for i in self.q_atom_ids:
            # if an atom ID is not part of the query clusters and not a descriptor
            if self.q_atom_ids[i] not in query_clusters and i not in self.q_atomg_nodeind_desc:
                if self.atomistic_dfs.node_info_1[0][i] == "":
                    continue
                self.query_endgroups.add(self.q_atom_ids[i]) 
                
        target_clusters = set()
        for i in self.target_clusters:
            target_clusters.update(i) 
        self.target_endgroups = set()
        for i in self.t_atom_ids:
            if self.t_atom_ids[i] not in target_clusters and i not in self.t_atomg_nodeind_desc:
                self.target_endgroups.add(self.t_atom_ids[i])     
        self.atomistic_dfs.set_clusters(self.query_clusters, self.target_clusters, self.q_macrocycle, self.t_macrocycle) 

    def search_repeats_endgroups(self): 
        # if there are no other objects, then the # of query objects = # target objects
        if self.no_other_objects:
            if len(self.query_cycles) != len(self.target_cycles):
                return False

        # the query must have <= number of target objects
        if len(self.query_cycles) > len(self.target_cycles): 
            return False
        
        # SMARTS in SMILES search, no stochastic objects
        if len(self.query_cycles) == len(self.target_cycles) == 0: 
            # a small molecule query cannot have more atoms than a small molecule target
            if self.q_atom_graph.number_of_nodes() > self.t_atom_graph.number_of_nodes():
                return False
            # initializes DFS, the two input parameters are which topological IDs in the query
            # are searched in which topological IDs in the target, for instance a repeat unit searched in a target stochastic object
            # for search, which for a small molecule would only be {1}
            self.atomistic_dfs.init_dfs({1},{1})
            for match in self.atomistic_dfs.generate_matches():
                # the very first match is True
                return True
            # if no match at at all, then False
            return False

        # SMARTS len(self.query_cycles) == 0 in BigSMILES search
        if len(self.query_cycles) == 0 and len(self.target_cycles) > 0: 
            # now the entire query {1} is searched in all topological IDs in the target
            self.atomistic_dfs.init_dfs({1}, self.all_ids_2)
            for match in self.atomistic_dfs.generate_matches():
                # print("match: ", match)
                if match != dict():
                    return True
            # if no match at at all, then False
            return False
        
        # BigSMARTS in BigSMILES
        # PARSE stochastic objects and non-atom elements, which we call local elements
        # ru_local would be elements like no other repeat units !*, no other end groups !*
        # for example, {[][<]CCO[>],!*[]}
        # {'wildcard_cluster': False, 'ru_local_el': {'!*'}, 'endgrp_local_el': set()}}
        # {[][<]CCO[>],!*;!*[]}
        # {'wildcard_cluster': False, 'ru_local_el': {'!*'}, 'endgrp_local_el': {'!*'}}}
        # {[][]}
        # {'wildcard_cluster': True, 'ru_local_el': set(), 'endgrp_local_el': set()}
        local_el_atom = nx.get_node_attributes(self.q_atom_graph, "local_el")
        self.local_el = dict()
        for atom in local_el_atom:
            cluster = [cl for cl in self.query_clusters if self.q_atom_ids[atom] in cl]
            if len(cluster) == 0:
                continue
            cluster = sorted(cluster[0])
            if tuple(cluster) in self.local_el:
                self.local_el[tuple(cluster)]["ru_local_el"].update(local_el_atom[atom]["ru_local_el"])
                self.local_el[tuple(cluster)]["endgrp_local_el"].update(local_el_atom[atom]["endgrp_local_el"])
            else:
                self.local_el[tuple(cluster)] = local_el_atom[atom]

        # First, search cycles or repeat units
        matches = dict() 
        # Iterate through each stochastic object, which I call query_cluster
        for query_cluster in self.query_cycles:        
            def search_smarts_in_clusters(smarts_list, target_cluster):
                # searches smarts in clusters
                for local in smarts_list:
                    if local == "!*":
                        continue
                    local_graphs = search_tools.generate_NetworkX_graphs(input_string = local, is_bigsmarts = True)  
                    target_graphs = search_tools.generate_NetworkX_graphs(input_string = self.bigsmiles, is_bigsmarts = False) 
                    local_search = atomistic.Atom_Graph_Matcher(local_graphs["atomistic"], target_graphs["atomistic"], local_graphs["descriptors"], target_graphs["descriptors"])
                    local_search.init_dfs({1}, target_cluster)
                    found_in_cluster = False
                    for match in local_search.generate_matches():
                        found_in_cluster = True
                        break
                    if not found_in_cluster:
                        return False
                return True

            def match_cluster_to_cluster(query_cluster, target_cluster):
                # match every atom in query cluster to target cluster
                hits = dict()
                for atom_q in self.q_atom_graph:
                    if self.q_atom_ids[atom_q] in query_cluster and atom_q not in self.q_atomg_nodeind_desc:        
                        atom_list = []                        
                        for atom_t in self.t_atom_graph:
                            if self.t_atom_ids[atom_t] in target_cluster and atom_t not in self.t_atomg_nodeind_desc:
                                atom_list.append(atom_t)
                        hits[atom_q] = atom_list
                return hits
            
            def is_wildcard_cycle(cycle):
                symbols = nx.get_node_attributes(self.q_atom_graph, "symbol")
                ids = nx.get_node_attributes(self.q_atom_graph, "ids")
                for key in ids:
                    if ids[key] in cycle and key not in self.q_atomg_nodeind_desc and symbols[key] not in ["*", "?*"]:
                        return False
                return True

            def match_cycle_to_cycle(query_cycle, target_cycle):
                # map every atom in query cycle to all atoms in target cycle
                hits = dict()
                for atom_q in self.q_atom_graph:
                    
                    # get all atoms in the query wildcard cycle
                    if self.q_atom_ids[atom_q] in query_cycle and atom_q not in self.q_atomg_nodeind_desc:
                        
                        # get all atoms in the target cycle
                        atom_list = []                        
                        for atom_t in self.t_atom_graph:
                            if self.t_atom_ids[atom_t] in target_cycle and atom_t not in self.t_atomg_nodeind_desc:
                                atom_list.append(atom_t)
                        hits[atom_q] = atom_list
                return hits
            
            cluster = set()
            for cyc in self.query_clusters:
                if query_cluster[0][0] in cyc:
                    cluster = sorted(cyc)

            # if the stochastic object is a wildcard object 
            if self.local_el[tuple(cluster)]["wildcard_cluster"]:
                hits = []
                # iterate through each stochastic object in the target
                for target_cluster in self.target_clusters:
                    # is a SMARTS is a localization element, it searches the target cluster
                    # {[]CCOCCO[]} there is a SMARTS in the stochastic object
                    if search_smarts_in_clusters(self.local_el[tuple(cluster)]["ru_local_el"], target_cluster):
                        # this matches a set of atoms in a query wildcard cluster to a set of atoms in a target cluster
                        hits.append(match_cluster_to_cluster(cluster, target_cluster))
                # map a cluster to a set of hits
                matches[tuple(cluster)] = hits
            else:    
                # iterate through each cycle of the query cluster
                for cycle in query_cluster: 
                    hits = []
                    if is_wildcard_cycle(cycle):
                        for target_cluster in self.target_cycles:
                            if search_smarts_in_clusters(self.local_el[tuple(cluster)]["ru_local_el"], target_cluster):
                                for target_cycle in target_cluster:
                                    hits.append(match_cycle_to_cycle(cycle, target_cycle))
                    else:
                        # for target_cluster in self.target_cycles:
                        #     for target_cycle in target_cluster:
                        #         self.atomistic_dfs.init_dfs(cycle, target_cycle)
                        #         found = False  
                        #         for match in self.atomistic_dfs.generate_matches():
                        #             values = list(match.values())[0]
                        #             if type(values) == list:
                        #                 values = values[0]
                        #             for target_cluster in self.target_clusters:
                        #                 if self.t_atom_ids[values] in target_cluster:
                        #                     if search_smarts_in_clusters(self.local_el[tuple(cluster)]["ru_local_el"], target_cluster):
                        #                         hits.append(match) 
                        #                         found = True
                        #                         break
                        #             if found:
                        #               break

                        # search cycle in all target topological IDs
                        self.atomistic_dfs.init_dfs(cycle, self.all_ids_2) 
                        # for each match
                        for match in self.atomistic_dfs.generate_matches():
                            values = list(match.values())[0]
                            if type(values) == list:
                                values = values[0]
                            # identify the stochastic object the query cycle matches to and search SMARTS
                            # {[][<]CCO[>],CCOCCO[]}
                            for target_cluster in self.target_clusters:
                                if self.t_atom_ids[values] in target_cluster:
                                    if search_smarts_in_clusters(self.local_el[tuple(cluster)]["ru_local_el"], target_cluster):
                                        hits.append(match) 

                    if hits == []: 
                        return False
                    # map a cycle to a list of hits
                    matches[tuple(cycle)] = hits 

        # Endgroup search   
        for id in self.query_endgroups:
            def is_wildcard_endgroup(endgroup):
                symbols = nx.get_node_attributes(self.q_atom_graph, "symbol")
                ids = nx.get_node_attributes(self.q_atom_graph, "ids") 
                atoms = []
                for key in ids:
                    if ids[key] == endgroup and key not in self.q_atomg_nodeind_desc and symbols[key] not in ["?*", "*"]:
                        atoms.append(symbols[key]) 
                if atoms == []:
                    return True
                return False

            def generate_endgroup_candidates(endgroup_query):
                target_ends = list(self.target_endgroups)
                # print("target_ends: ", target_ends)

                # get number of adjacent descriptors in query endgroup
                ids_topology_q = nx.get_node_attributes(self.q_top_undir_graph, "ids")
                num_desc = 0
                for key in ids_topology_q:
                    if ids_topology_q[key] == endgroup_query:
                        ids_desc = []
                        for i in self.q_atomg_nodeind_desc:
                            ids_desc.append(ids_topology_q[i])
                        for neighbor in self.q_top_undir_graph[key]:
                            if ids_topology_q[neighbor] in ids_desc:
                                num_desc += 1

                # print("num desc: ", num_desc, target_ends) 
                # suggest target endgroups that have at least that many descriptors
                endgroup_combos = []
                for i in range(len(target_ends) + 1):
                    for j in itertools.combinations(range(len(target_ends)), i):
                        j = list(j)
                        ids_topology_t = nx.get_node_attributes(self.t_top_undir_graph, "ids")
                        ids_desc = [] 
                        for k in self.t_atomg_nodeind_desc:
                            ids_desc.append(ids_topology_t[k])
                        # print("ids_desc: ", ids_desc)
                        ids_eg = [] 
                        for k in j:
                            ids_eg.append(target_ends[k])

                        def valid_endgroup(j):
                            for m in range(len(j)):
                                for n in range(m + 1, len(j)):
                                    def is_valid(path):
                                        for p in path:
                                            if ids_topology_t[p] not in ids_desc and ids_topology_t[p] not in ids_eg:
                                                return False
                                        return True
                                    
                                    m_node_id = -1
                                    n_node_id = -1
                                    for key in ids_topology_t:
                                        if ids_topology_t[key] == target_ends[j[m]]:
                                            m_node_id = key
                                        if ids_topology_t[key] == target_ends[j[n]]:
                                            n_node_id = key
                                
                                    valid = False
                                    for path in nx.all_simple_paths(self.t_top_undir_graph, m_node_id, n_node_id):
                                        if is_valid(path):
                                            valid = True
                                            break
                                    if not valid:
                                        return False

                            return True

                        if valid_endgroup(j):
                            endgroup_combos.append([target_ends[k] for k in j])

                endgroup_combos.remove([])  
                # print("endgroup combos: ", endgroup_combos)
                endgroups = []
                for i in range(len(endgroup_combos)):          
                    def count_desc(endgroup_combo):
                        ids_topology_t = nx.get_node_attributes(self.t_top_undir_graph, "ids")
                        num_desc_t = set()
                        for key in ids_topology_t:
                            if ids_topology_t[key] in endgroup_combo:
                                ids_desc = []
                                for j in self.t_atomg_nodeind_desc:
                                    ids_desc.append(ids_topology_t[j])
                                
                                for neighbor in self.t_top_undir_graph[key]:
                                    if ids_topology_t[neighbor] in ids_desc:
                                        num_desc_t.add(neighbor)

                        return len(num_desc_t)

                    if count_desc(endgroup_combos[i]) >= num_desc:
                        endgroups.append(endgroup_combos[i])

                return endgroups

            def match_endgroup_to_endgroup(endgroup_query, endgroup_target):
                # get all atoms in query endgroup
                atom_list = []
                for atom in self.q_atom_graph:
                    if self.q_atom_ids[atom] == endgroup_query and atom not in self.q_atomg_nodeind_desc:
                        atom_list.append(atom)
                
                # when counting neighboring stochastic objects to endgroups in query and target, 
                # make sure there is only one connection to atoms
                count_q = 0
                for desc in self.q_atomg_nodeind_desc:
                    count = 0
                    for neighbor in self.q_atom_graph[desc]:
                        if neighbor in atom_list:
                            count += 1
                    if count == 1:
                        count_q += 1
                atom_list = []
                for atom in self.t_atom_graph:
                    if self.t_atom_ids[atom] in endgroup_target and atom not in self.t_atomg_nodeind_desc:
                        atom_list.append(atom) 
                
                count_t = 0
                for desc in self.t_atomg_nodeind_desc:
                    count = 0
                    for neighbor in self.t_atom_graph[desc]:
                        if neighbor in atom_list:
                            count += 1
                    if count == 1:
                        count_t += 1
                if count_q > count_t:
                    return []

                for i in endgroup_target:
                    ids_undir = nx.get_node_attributes(self.t_top_undir_graph, "ids")
                    inv_map = {v: k for k, v in ids_undir.items()}
                    if len(self.t_top_undir_graph[inv_map[i]]) == 1:
                        return []

                # map every atom in endgroup_query to all endroups in endgroup
                hits = dict()
                for atom in self.q_atom_graph:
                    if self.q_atom_ids[atom] == endgroup_query and atom not in self.q_atomg_nodeind_desc:
                        atom_list = []                        
                        for atom_t in self.t_atom_graph:
                            if self.t_atom_ids[atom_t] in endgroup_target and atom_t not in self.t_atomg_nodeind_desc:
                                atom_list.append(atom_t)
                        hits[atom] = atom_list 

                return hits
                
            hits = []
            if is_wildcard_endgroup(id):
                endgroups = generate_endgroup_candidates(id) 
                for e in endgroups:     
                    h = match_endgroup_to_endgroup(id, e)
                    if h != []:
                        hits.append(h)
            else:
                self.atomistic_dfs.init_dfs({id}, self.all_ids_2)
                for match in self.atomistic_dfs.generate_matches():
                    hits.append(match)
                if hits == []:
                    return False
            
            # map an end group topological ID to a list of hits
            matches[(id,)] = hits
        
        # store a list of matches
        self.ru_eg_matches = matches
        self.init_top_dfs()  
        # traverse through 
        for match in self.level():
            return True
        return False
            
    def init_top_dfs(self): 
        # atom mapping of key --> value
        self.ovrll_atm_mapped = {}

        # atom mapping of key --> depth
        self.atm_mapped_order = {}  

        # list of smiles that have already been mapped
        self.smiles_mapped = set()

        # next set of smiles to choose from
        self.smiles_next = {} 

        # state of mapping
        self.state = GMState(self) 

    def level(self):
        if len(self.smiles_mapped) == len(self.ru_eg_matches):
            def checks():
                # are wildcard cycles adjacent
                def wildcard_adjacency():
                    for key1 in self.ovrll_atm_mapped:
                        if self.q_atom_ids[key1] in self.query_endgroups and self.atomistic_dfs.node_info_1[0][key1] == "?*":
                            adjacent = []
                            for key2 in self.ovrll_atm_mapped:
                                if descriptors_separating(key1, key2, self.atomistic_dfs.bonds_1, self.atomistic_dfs.descriptors_1):
                                    adjacent.append(key2)
                            for atom in adjacent:
                                if type(self.ovrll_atm_mapped[atom]) == list:
                                    reference = self.t_atom_ids[self.ovrll_atm_mapped[atom][0]]
                                else:
                                    reference = self.t_atom_ids[self.ovrll_atm_mapped[atom]]
                                for t in self.target_clusters:
                                    if reference in t:
                                        reference_cluster = t
                                found = False
                                for key3 in self.ovrll_atm_mapped[key1]:
                                    for key4 in self.t_atom_ids:
                                        if self.t_atom_ids[key4] in reference_cluster:
                                            if descriptors_separating(key3, key4, self.atomistic_dfs.bonds_2, self.atomistic_dfs.descriptors_2):
                                                found = True
                                if not found:
                                    return False
                    for cluster in self.query_cycles:
                        for cycle in cluster:
                            ids_q = nx.get_node_attributes(self.q_atom_graph, "ids")
                            ids_t = nx.get_node_attributes(self.t_atom_graph, "ids")
                            count_wildcard = 0
                            count_desc = 0
                            for key in ids_q:
                                if ids_q[key] in cycle:
                                    if self.atomistic_dfs.node_info_1[0][key] == "?*":
                                        count_wildcard += 1
                                    if key in self.atomistic_dfs.descriptors_1:
                                        count_desc += 1
                            seg = count_wildcard == 2 and count_desc == 2
                            graft = count_wildcard == 3 and count_desc == 1
                            if seg or graft:
                                nested_cycle = []
                                nonnested_cycle = []
                                level = nx.get_node_attributes(self.q_atom_graph, "level")
                                base = cycle
                                for key in level:
                                    if level[key] == 1:
                                        for other_cycle in cluster:
                                            if ids_q[key] in other_cycle:
                                                nested = other_cycle
                                                break
                                level = nx.get_node_attributes(self.t_atom_graph, "level")
                                for key in self.ovrll_atm_mapped:
                                    if ids_q[key] in nested and key not in self.q_atomg_nodeind_desc:
                                        mapped = self.ovrll_atm_mapped[key]
                                        if type(mapped) == list:
                                            for m in mapped:
                                                if m not in self.t_atomg_nodeind_desc:
                                                    chosen = m
                                        else:
                                            chosen = mapped
                                        for cluster2 in self.target_cycles:
                                            for cycle2 in cluster2:
                                                if ids_t[chosen] in cycle2 and level[chosen] == 1:
                                                    nested_cycle = cycle2
                                                    break
                                    if ids_q[key] in base and key not in self.q_atomg_nodeind_desc:
                                        mapped = self.ovrll_atm_mapped[key]
                                        if type(mapped) == list:
                                            for m in mapped:
                                                if m not in self.t_atomg_nodeind_desc:
                                                    chosen = m
                                        else:
                                            chosen = mapped
                                        for cluster2 in self.target_cycles:
                                            for cycle2 in cluster2:
                                                if ids_t[chosen] in cycle2 and level[chosen] == 0:
                                                    nonnested_cycle = cycle2
                                                    break
                                shared = any(i in nonnested_cycle for i in nested_cycle)                         
                                if seg and shared:
                                    return True
                                if seg and not shared:
                                    return False
                                if graft and shared:
                                    return False
                                if graft and not shared:
                                    bonds = nx.get_edge_attributes(self.t_top_undir_graph, "bond_type")
                                    ids = nx.get_node_attributes(self.t_top_undir_graph, "ids")
                                    for g in bonds:
                                        ids1 = ids[g[0]]
                                        ids2 = ids[g[1]]
                                        if ids1 in nonnested_cycle and ids2 in nested_cycle or \
                                            ids2 in nonnested_cycle and ids1 in nested_cycle:
                                            return True
                                    return False
                    return True
                # do nested query cycles match to nested target cycles
                def nested_nested():
                    level_1 = nx.get_node_attributes(self.q_atom_graph, "level")
                    level_2 = nx.get_node_attributes(self.t_atom_graph, "level")
                    for key in level_1:
                        if level_1[key] == 1:
                            mapped = self.ovrll_atm_mapped[key]
                            if type(mapped) == list:
                                mapped = mapped[0]
                            if level_2[mapped] == 0:
                                return False
                    return True
                
                def cluster_cluster():
                    target_cluster = []
                    for q_cluster in self.query_clusters:

                        # get all query atoms from the same cluster
                        query_atoms = []
                        for atom in self.ovrll_atm_mapped:
                            if self.q_atom_ids[atom] in q_cluster:
                                query_atoms.append(atom)
                                            
                        for i in range(len(self.target_clusters)):
                            t_cluster = self.target_clusters[i]
                            
                            # all query_atoms in that cluster are found in a single target cluster
                            count = 0
                            for atom in query_atoms:
                                if type(self.ovrll_atm_mapped[atom]) == list:
                                    a = self.ovrll_atm_mapped[atom][0]
                                    if self.t_atom_ids[a] in t_cluster:
                                        count += 1
                                else:
                                    if self.t_atom_ids[self.ovrll_atm_mapped[atom]] in t_cluster:
                                        count += 1

                            if count == len(query_atoms):
                                if t_cluster in target_cluster:
                                    return False
                                else:
                                    target_cluster.append(t_cluster)
                                break
                            elif i == len(self.target_clusters) - 1:
                                return False
                        
                    return True

                def endgroup_matches():
                    # does the same target endgroup atom map more than once?
                    for key in self.t_atom_ids:
                        if self.t_atom_ids[key] in self.target_endgroups:
                            count = 0
                            for key2 in self.ovrll_atm_mapped:
                                if type(self.ovrll_atm_mapped[key2]) == int and key == self.ovrll_atm_mapped[key2]:
                                    count += 1
                                if count > 1:
                                    return False
                    already_mapped = []
                    for eg in self.query_endgroups:
                        for key in self.ovrll_atm_mapped:
                            if type(self.ovrll_atm_mapped[key]) == list and self.q_atom_ids[key] == eg:
                                if self.ovrll_atm_mapped[key] in already_mapped:
                                    return False
                                already_mapped.append(self.ovrll_atm_mapped[key])
                                break
                    return True

                def exclamation_element():
                    # !* in repeat units
                    #   cycle to cycle mapping, for each query cycle, no target node is repeated
                    #   all ids in the target cluster must be hit
                    # !* in endgroups
                    #   all endgroup and cycle ids outside the target cluster must be hit

                    # get set of target cycles
                    target_cycles = []
                    for t_atom in self.target_cycles:
                        for b in t_atom:
                            target_cycles.append(b)

                    # iterate through each query cluster and determine if "!*" is present
                    for q_cluster in self.local_el:

                        if "!*" in self.local_el[q_cluster]["ru_local_el"]:     
                            # store the target cluster that the query cluster maps to
                            cycles_mapped = dict()
                            target_cycles_mapped = []

                            for q_cycle_set in self.query_cycles:
                                for q_cycle in q_cycle_set:
                                    if q_cycle[0] in q_cluster:

                                        # for each cycle, get all query atoms and atoms mapped in that cycle
                                        atoms = []
                                        atoms_mapped = []
                                        for t_atom in self.ovrll_atm_mapped:
                                            if self.q_atom_ids[t_atom] in q_cycle:
                                                atoms.append(t_atom)
                                                if type(self.ovrll_atm_mapped[t_atom]) == list:
                                                    atoms_mapped.append(tuple(self.ovrll_atm_mapped[t_atom]))
                                                else:
                                                    atoms_mapped.append(self.ovrll_atm_mapped[t_atom])

                                        # no node in atoms_mapped can be repeated
                                        if len(atoms_mapped) != len(set(atoms_mapped)):
                                            return False

                                        # only 1 cycle in the target allowed to be mapped
                                        mapped = []
                                        for t_atom in atoms_mapped:
                                            for t_cycle in target_cycles:
                                                if self.t_atom_ids[t_atom] in t_cycle:
                                                    if t_cycle != mapped and mapped != []:
                                                        return False
                                                    if t_cycle[0] not in self.target_endgroups:
                                                        mapped = t_cycle
                                                    break  

                                        cycles_mapped[tuple(q_cycle)] = mapped
                                        target_cycles_mapped.append(sorted(mapped))
                                        
                            # all cycles in the target cluster should be mapped
                            for cluster in self.target_cycles:         
                                # iterate through all target cycles in each cluster
                                # check if every cycle in target cycles is found in target cycles mapped
                                t_cycles = [sorted(t_cycle) for t_cycle in cluster]
                                if target_cycles_mapped[0] in t_cycles:
                                    for r in t_cycles:
                                        if r not in target_cycles_mapped:
                                            return False

                        if "!*" in self.local_el[q_cluster]["endgrp_local_el"]:
                            for key in self.ovrll_atm_mapped:
                                # iterate through every query atom in cluster and get mapped atom
                                if self.q_atom_ids[key] in q_cluster:
                                    mapped = self.ovrll_atm_mapped[key]
                                    if type(mapped) == list:
                                        mapped = mapped[0]
                                    
                                    # find target cluster than query cluster maps to
                                    for cluster in self.target_clusters:
                                        if self.t_atom_ids[mapped] in cluster:
                                            # if there is any atom not part of the cluster, return False
                                            for key2 in self.t_atom_ids:
                                                if self.t_atom_ids[key2] not in cluster:
                                                    return False
                                    break

                    return True

                return endgroup_matches() and cluster_cluster() and wildcard_adjacency() and nested_nested() and exclamation_element()
            
            if checks():
                # print("Final: ", self.ovrll_atm_mapped)
                yield self.ovrll_atm_mapped 
        
        else: 
            for smiles_chosen, mapping_suggested in self.suggest_mappings():
                if self.are_mappings_adjacent(smiles_chosen, mapping_suggested):
                    newstate = self.state.__class__(self, smiles_chosen, mapping_suggested)
                    yield from self.level() 
                    newstate.backtrack()  

    def suggest_mappings(self):
        T1_inout = [node for node in self.smiles_next if node not in self.smiles_mapped]
        if T1_inout: 
            smiles = T1_inout[0]
        else:
            smiles = tuple(self.ru_eg_matches.keys())[0]

        mappings = self.ru_eg_matches[smiles]
        for mapping_suggested in mappings:
            yield smiles, mapping_suggested

    def are_mappings_adjacent(self, smiles_chosen, mapping_suggested):
        current_cluster = set()
        for cyc_endgroup in self.smiles_mapped:
            for id in cyc_endgroup:
                current_cluster.add(id)

        for key in smiles_chosen:
            current_cluster.add(key)
        
        mapping = copy.deepcopy(self.ovrll_atm_mapped)
        for key in mapping_suggested:
            mapping[key] = mapping_suggested[key]
        
        return enforce_adjacency(mapping, self.atomistic_dfs, current_cluster)

class GMState:   
    def __init__(self, GM, smiles_chosen=None, mapping_suggested=None):
        self.GM = GM

        # Initialize the last stored node pair.
        self.smiles_node = None
        self.depth = len(GM.smiles_next) 

        if smiles_chosen is None:
            GM.ovrll_atm_mapped = {}
            GM.atm_mapped_order = {}
            GM.smiles_mapped = set()
            GM.smiles_next = {}
    
        if smiles_chosen is not None:
            # keep track of the smiles that was chosen
            self.smiles_node = smiles_chosen 
            self.depth = len(GM.smiles_mapped) 

            # add mapping suggested to overall mapping
            for atom in mapping_suggested:
                if atom not in GM.ovrll_atm_mapped:
                    GM.ovrll_atm_mapped[atom] = mapping_suggested[atom]
                    GM.atm_mapped_order[atom] = self.depth
            
            # add smiles that was just chosen to smiles mapped
            GM.smiles_mapped.add(smiles_chosen)

            # augment smiles_next
            if smiles_chosen not in GM.smiles_next:
                GM.smiles_next[smiles_chosen] = self.depth
            
            # Get repeat unit keys to determine next cycle or endgroup to search
            ru_keys = list(GM.ru_eg_matches.keys())

            # Get topological undirected graph representation to run DFS on topological graph
            q_top_undir_graph_ids = nx.get_node_attributes(GM.q_top_undir_graph, "ids")

            # iterate through all ids in mapped smiles
            for atom in GM.smiles_mapped.copy():
                for id in atom:

                    # get node in contracted graph
                    for node in q_top_undir_graph_ids:
                        if id == q_top_undir_graph_ids[node]:

                            # get neighbors of the node in contracted graph
                            for neighbor in GM.q_top_undir_graph[node]:

                                # there must be an adjacent cycle or endgroup
                                found_adj_cyc_endgrp = False 

                                # get all neighbors for that node
                                for k in ru_keys: 
                                    if q_top_undir_graph_ids[neighbor] in k:
                                        
                                        found_adj_cyc_endgrp = True
                                        
                                        if tuple(k) not in GM.smiles_next:

                                            # get adjacent cycles and endgroups
                                            GM.smiles_next[tuple(k)] = self.depth
                                
                                # if a descriptor is landed on but there is no cycle or endgroup
                                if not found_adj_cyc_endgrp:
                                    for neighbor2 in GM.q_top_undir_graph[neighbor]:

                                        for k in ru_keys: 
                                            if q_top_undir_graph_ids[neighbor2] in k:
                                                
                                                if tuple(k) not in GM.smiles_next:

                                                    # get adjacent cycles and endgroups
                                                    GM.smiles_next[tuple(k)] = self.depth

    def backtrack(self):
        for key in self.GM.ovrll_atm_mapped.copy():
            if self.GM.atm_mapped_order[key] == self.depth:
                del self.GM.ovrll_atm_mapped[key]
                del self.GM.atm_mapped_order[key]

        if self.smiles_node is not None:
            self.GM.smiles_mapped.remove(self.smiles_node)

        for node in self.GM.smiles_next.copy():
            if self.GM.smiles_next[node] == self.depth:
                del self.GM.smiles_next[node] 