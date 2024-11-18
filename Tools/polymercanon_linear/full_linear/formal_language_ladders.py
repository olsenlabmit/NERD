import re
import rdkit
from rdkit import Chem

import formal_language_non_ladders 
import priority_rules_non_ladders   

def generate_transitions(bigsmiles):
    bigsmiles = bigsmiles[1:-1]
    d = re.findall(r"\[.\d+\[.\d+\]\d+\]", bigsmiles)
    
    descriptors = []
    alphabets = ["[*:1]","[*:2]","[*:3]","[*:4]"]
    for i in range(len(d)):
        outer = re.findall(r"\[.\d+", d[i])[0][1:]
        inner = re.findall(r"\[.\d+", d[i])[1][1:]
        group_id = re.findall(r"\d+\]", d[i])[1][:-1]
        bigsmiles = bigsmiles.replace(d[i],alphabets[i],1)
        descriptors.append([d[i], outer, inner, group_id])

    inner = [descriptors[0][2],
                descriptors[1][2],
                descriptors[2][2],
                descriptors[3][2]]
                
    horizontal = formal_language_non_ladders.get_compatible_descriptor(inner[0]) == inner[2] and \
                formal_language_non_ladders.get_compatible_descriptor(inner[1]) == inner[3]
    vertical = formal_language_non_ladders.get_compatible_descriptor(inner[0]) == inner[3] and \
                formal_language_non_ladders.get_compatible_descriptor(inner[1]) == inner[2]
    
    new = bigsmiles
    orig = bigsmiles
    if horizontal and vertical:
        # flip descriptors (1-->2 and 3-->4) to determine if repeat unit is symmetric across horizontal axis
        new = new.replace("[*:1]","X")
        new = new.replace("[*:2]","[*:1]") 
        new = new.replace("X","[*:2]")    
        new = new.replace("[*:3]","X")
        new = new.replace("[*:4]","[*:3]") 
        new = new.replace("X","[*:4]")  
    
        if Chem.MolToSmiles(Chem.MolFromSmiles(new)) == Chem.MolToSmiles(Chem.MolFromSmiles(orig)):
            bigsmiles = "{[][<]C[>][]}"
            bigsmiles = formal_language_non_ladders.process(bigsmiles)
            automata, backbone, o = formal_language_non_ladders.generate_paths(bigsmiles, "start", "end", [], dict(), 1)
            initialState, finalState, statesList, alphabet, transition, alphabet_dictionary = formal_language_non_ladders.generate_transitions(automata, backbone)
            alphabet_dictionary[new] = alphabet_dictionary.pop("C([*:1])[*:2]")
        else:
            bigsmiles = "{[][<]C[>],[<]O[>][]}"
            bigsmiles = formal_language_non_ladders.process(bigsmiles)
            automata, backbone, o = formal_language_non_ladders.generate_paths(bigsmiles, "start", "end", [], dict(), 1)
            initialState, finalState, statesList, alphabet, transition, alphabet_dictionary = formal_language_non_ladders.generate_transitions(automata, backbone)
            alphabet_dictionary[new] = alphabet_dictionary.pop("C([*:1])[*:2]")
            alphabet_dictionary[orig] = alphabet_dictionary.pop("O([*:1])[*:2]")
    
    return initialState, finalState, statesList, alphabet, transition, alphabet_dictionary