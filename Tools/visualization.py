from rdkit import Chem
from rdkit.Chem import Draw
import networkx as nx
from pygraphviz import *

def visualize_alphabets(alphabet_dictionary):
    mols = []
    legend = []
    for key in alphabet_dictionary:
        mols.append(Chem.MolFromSmiles(key))
        legend.append(alphabet_dictionary[key])
    img = Draw.MolsToGridImage(mols, molsPerRow=6, subImgSize=(200, 250), useSVG=True, legends=legend)
    filename = "graph_outputs/alphabets.svg"
    with open(filename, 'w') as f:
        f.write(img)   

def visualize_bu(polymer, name):
    isStart = nx.get_node_attributes(polymer, "isStart")
    isEnd = nx.get_node_attributes(polymer, "isEnd")
    isBranch = nx.get_node_attributes(polymer, "isBranch")
    state = nx.get_node_attributes(polymer, "state")
    depiction_id = dict()
    add_on = ""
    count = 1
    for key in isBranch:
        add_on += str(key) + "[shape = plaintext, color = white, style = filled];"
    for key in state:
        if isStart[key]:
            add_on += "I" + str(count) + "[shape=point, color=white];"
            add_on += "I" + str(count) + "->" + str(key) + "[color=darkgreen, penwidth=4];"
            count += 1
        if isEnd[key]:
            add_on += str(key) + "[shape = doublecircle, color = red, style = filled];"
        if not isStart[key] and not isEnd[key]:
            depiction_id[key] = " "
    nx.set_node_attributes(polymer, depiction_id, "depiction_id")
    
    visual = nx.nx_agraph.write_dot(polymer, "graph_outputs/dot/" + name + "_dot.txt")
    text_file = open("graph_outputs/dot/" + name + "_dot.txt", "r")
    data = text_file.read()
    data = data.replace("digraph \"\" {", "digraph \"\" {" + add_on)
    data = data.replace("alpha", "label")
    data = data.replace("depiction_id=", "shape=circle, label=")
    text_file.close() 
    G = AGraph(data) 
    G.draw("graph_outputs/" + name + ".svg", prog='dot')
    return polymer

def set_colors(graphs, graph_name):
    colorscheme = "set35"
    colors = [5,4,3,2,1]
    style = "filled"
    
    symbols = nx.get_node_attributes(graphs[graph_name], "symbol")
    ids = nx.get_node_attributes(graphs[graph_name], "ids")
    groups = {}
    for key in ids:
        if ">" in symbols[key] or "<" in symbols[key] or "$" in symbols[key]:
            continue
        if ids[key] in groups:
            groups[ids[key]].append(key)
        else:
            groups[ids[key]] = [key]
    groups = list(groups.values())
    counter = 0
    legend = {}
    for g in groups:
        legend[tuple(g)] = {"colorscheme": colorscheme, "color": colors[counter], "style": style}
        counter += 1
        counter = counter % len(colors)
    
    color = {}
    colorscheme = {}
    style = {}
    for key in symbols:
        for key2 in legend:
            if key in key2:
                color[key] = legend[key2]["color"]
                colorscheme[key] = legend[key2]["colorscheme"]
                style[key] = legend[key2]["style"]
    nx.set_node_attributes(graphs[graph_name], color, "color") 
    nx.set_node_attributes(graphs[graph_name], colorscheme, "colorscheme") 
    nx.set_node_attributes(graphs[graph_name], style, "style") 
    return graphs

def visualize_NetworkX_graphs(graphs, graph_name, file_name):
    symbols = nx.get_node_attributes(graphs["atomistic"], "symbol")
    level = nx.get_node_attributes(graphs["atomistic"], "level")
    node_id = dict()
    for key in symbols:
        # node_id[key] = symbols[key]
        node_id[key] = str(key) + " " + symbols[key] # + " " + str(graphs["ids"][key])
        # if key in level:
        #     node_id[key] = str(key) + " " + str(level[key]) + " " + symbols[key]
        # else:
        #     node_id[key] = str(key) + " " + symbols[key]
    nx.set_node_attributes(graphs["atomistic"], node_id, "depiction_id")

    symbols = nx.get_node_attributes(graphs["top_undir"], "symbol")
    node_id = dict()
    for key in graphs["ids"]:
        if key in graphs["descriptors"]:
            node_id[key] = symbols[key] + " " + str(graphs["ids"][key])
        else:
            node_id[key] = str(graphs["ids"][key]) 
    nx.set_node_attributes(graphs["top_undir"], node_id, "depiction_id")

    node_id = dict()
    for key in graphs["ids"]:
        if key in graphs["descriptors"]:
            node_id[key] = symbols[key] + " " + str(graphs["ids"][key])
            # node_id[key] = str(key) + " " + symbols[key] + " " + str(graphs["ids"][key])
        else:
            node_id[key] = str(graphs["ids"][key])
            # node_id[key] = str(key) + " " + str(graphs["ids"][key])
    nx.set_node_attributes(graphs["topology"], node_id, "depiction_id")

    node_id = dict()
    for key in graphs["ids"]:
        if key in graphs["descriptors"]:
            node_id[key] = symbols[key] + " " + str(graphs["ids"][key])
        else:
            node_id[key] = str(graphs["ids"][key])
    nx.set_node_attributes(graphs["multidigraph"], node_id, "depiction_id")

    visual = nx.nx_agraph.write_dot(graphs[graph_name], "graph_outputs/dot/" + file_name + "_dot.txt")
    text_file = open("graph_outputs/dot/" + file_name + "_dot.txt", "r")
    data = text_file.read()
    data = data.replace("depiction_id=", "shape=circle, fontsize=\"30pt\", label=")
    if graph_name == "atomistic": 
        data = data.replace("bond_type", "fontsize=\"25pt\", label")
    text_file.close() 
    G = AGraph(data) 
    G.draw("graph_outputs/" + file_name + ".svg", prog='dot')