import copy
import networkx as nx

def convert_tuple_string(tup):
    st = ','.join(map(str, tup))
    return st

def get_comp_descriptor(desc):
    if "1" in desc:
        return "2" + desc[1:]
    return "1" + desc[1:]

def get_node_index(node):
    return int(node[:node.index("_")])

def get_iteration(node):
    return node[node.index("_") + 1:]

def feed_to_bonds_n(one, two):
    return tuple(sorted([one, two]))

def directly_connected(i, j, G):
    return i in G[j]

def descriptors_separating(atom_1, atom_2, bonds, descriptors):        
    d = []
    for desc in descriptors:
        try:
            pair_1 = feed_to_bonds_n(atom_1, desc)
            pair_2 = feed_to_bonds_n(atom_2, desc)
            a = "1" in bonds[pair_1] and "2" in bonds[pair_2]
            b = "2" in bonds[pair_1] and "1" in bonds[pair_2]
            if a or b:
                d.append(desc)
        except:
            continue
    return d

def enforce_adjacency(mapping, GM, q_cluster):  
    descriptors = []
    for id in q_cluster:
        # identify the descriptors in the cluster
        for key, value in GM.ids_1.items():
            if value == id:
                if key in GM.descriptors_1:
                    descriptors.append(key)
    
    for descriptor in descriptors:
        neighbors = GM.G1[descriptor]
        
        for node in neighbors:
            if GM.node_info_1[0][node] == "?*": 
                continue
            if GM.ids_1[node] in q_cluster:
                # [list of covalent connections to node], [list of connections across descriptor]
                direct, across = direct_across_connections(node, [descriptor, q_cluster], GM)
                for key in across:
                    across = across[key]

                node_mapped = mapping[node]                
                direct_mapped = [mapping[i] for i in direct]
                across_mapped = [mapping[i] for i in across]
                # [list of covalent connections to node mapped], {desc: list of connections across that descriptor}
                direct_available, across_available = direct_across_connections(node_mapped, [], GM)   
                
                for single in across_mapped:
                    # pick a single atom or path across the descriptor
                    mapped_1 = direct_mapped + [single]
                    mapped_2 = [copy.deepcopy(direct_available), copy.deepcopy(across_available)]
                    if not connections_found(mapped_1, mapped_2):
                        return False

    return True

def direct_across_connections(node, query, graphmatcher):

    if type(node) == list:
        return node, node
    if query != []:
        descriptor = query[0]
        cluster = query[1]

        descriptors = graphmatcher.descriptors_1
        graph = graphmatcher.G1
        ids = graphmatcher.ids_1
        bonds = graphmatcher.bonds_1
        all_clusters = graphmatcher.query_clusters

        empty_nodes = []
        symbols = nx.get_node_attributes(graph, "symbol")
        for atom in symbols:
            if symbols[atom] == "":
                empty_nodes.append(atom)
    else:
        descriptor = None
        cluster = set(graphmatcher.ids_2.values())

        descriptors = graphmatcher.descriptors_2
        graph = graphmatcher.G2
        ids = graphmatcher.ids_2
        bonds = graphmatcher.bonds_2
        all_clusters = graphmatcher.target_clusters

        empty_nodes = []
        symbols = nx.get_node_attributes(graph, "symbol")
        for atom in symbols:
            if symbols[atom] == "":
                empty_nodes.append(atom)

    direct = []
    across = dict()
    # iterate through all neighbors of the node, add to direct and across
    node_neighbors = list(graph[node])

    for i in range(len(node_neighbors)):
        if node_neighbors[i] not in descriptors:
            direct.append(node_neighbors[i])
        else: 
            if node_neighbors[i] == descriptor or descriptor is None and node_neighbors[i] in descriptors and ids[node_neighbors[i]] in cluster:
                across[node_neighbors[i]] = []
                pair_1 = feed_to_bonds_n(node, node_neighbors[i])
                for j in graph[node_neighbors[i]]:
                    if j in empty_nodes:
                        node_in = [c for c in all_clusters if ids[node] in c]
                        j_in = [c for c in all_clusters if ids[node_neighbors[i]] in c]
                        if node_in == j_in and node_in != []:
                            for k in graph[j]:
                                if k in descriptors and k != node_neighbors[i]:
                                    for atom in graph[k]:
                                        if atom in [descriptors + empty_nodes]:
                                            continue
                                        pair_2 = feed_to_bonds_n(k, atom)
                                        if bonds[pair_1] == get_comp_descriptor(bonds[pair_2]) and ids[atom] in cluster:
                                            across[node_neighbors[i]].append(atom) 
                    else:
                        pair_2 = feed_to_bonds_n(node_neighbors[i], j)
                        if bonds[pair_1] == get_comp_descriptor(bonds[pair_2]) and ids[j] in cluster:
                            across[node_neighbors[i]].append(j) 
    
    # if any element in direct or across has ?* in it, then delete it
    if query != []:
        direct = [atom for atom in direct if graphmatcher.node_info_1[0][atom] != "?*"]
        for descriptor in across:
            across[descriptor] = [atom for atom in across[descriptor] if graphmatcher.node_info_1[0][atom] != "?*"]
    return direct, across

def connections_found(mapped_1, mapped_2):
    num = dict()
    for node in mapped_1:
        if type(node) == list: 
            # if node is a path, only one atom in the list has to be equal to any atom in mapped_2 or in any list in mapped_2
            found = False
            for direct in mapped_2:
                for element in node:
                    if element in direct:
                        mapped_2.remove(direct)
                        found = True
            if found:
                continue

            for desc in mapped_2[1]:
                for element in node:
                    if element in direct:
                        found = True
                        mapped_2[1][desc].remove(node)
                        if desc not in num:
                            num[desc] = 1
                        else:
                            num[desc] += 1
                        if num[desc] > 1:
                            return False
            if not found:
                return False
        
        else:    
            if node in mapped_2[0]:
                mapped_2[0].remove(node) 
            else:
                found = False
                for desc in mapped_2[1]:
                    if node in mapped_2[1][desc]:
                        found = True
                        mapped_2[1][desc].remove(node)
                        if desc not in num:
                            num[desc] = 1
                        else:
                            num[desc] += 1
                        if num[desc] > 1:
                            return False
                if not found:
                    return False   

    return True