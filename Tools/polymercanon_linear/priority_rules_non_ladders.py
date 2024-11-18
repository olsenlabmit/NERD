import rdkit
from rdkit import Chem
import networkx as nx
import re, copy, collections
import pandas 

def process(bigsmiles):
    bigsmiles = bigsmiles.replace("}[H]","}")
    bigsmiles = bigsmiles.replace("[H]{","{")
    bigsmiles = bigsmiles.replace("[<]","[<1]")
    bigsmiles = bigsmiles.replace("[>]","[>1]")
    bigsmiles = bigsmiles.replace("[$]","[$1]")
    return bigsmiles

def get_objects(ends):
    i = 0
    o = []
    start = []
    while i < len(ends):
        if ends[i] == "{":
            start.append(i)
            object = "{"
            count = 1
            i += 1
            while count != 0:
                if ends[i] == "{":
                    count += 1
                if ends[i] == "}":
                    count -= 1
                object += ends[i]
                i += 1
            o.append(object)
        else:
            i += 1
    return [o, start]

def get_repeats(object):
    object = object[object.find("]")+1:]
    object = object[:object.rfind("[")]

    i = 0
    o = []
    count = 0
    repeat = ""
    object += ","
    while i < len(object):
        if object[i] == "{":
            count += 1
        if object[i] == "}":
            count -= 1
        if object[i] == "," and count == 0:
            o.append(repeat)
            repeat = ""
        else:
            repeat += object[i]
        i += 1
    return o

def sub_obj_Bk(bigsmiles):
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

def add_descriptors_back(d, index, canonical_order, smi):
    m = Chem.MolFromSmiles(smi)
    counter = 0
    iter = 0
    new_pos = []
    for i in index:
        new_pos.append(canonical_order.index(i))
    for q, atom in enumerate(m.GetAtoms()):
        smi_lowercase = smi.lower()
        # iterate through each atom until the a location is found of where the descriptors should be (a and b)
        while smi_lowercase[iter:iter + len(atom.GetSymbol())] != atom.GetSymbol().lower():
            iter += 1
        for i in range(len(new_pos)):
            if counter == new_pos[i]:
                smi = smi[0:iter - 1] + d[i] + smi[iter + 1 + len(atom.GetSymbol()):]
                iter += len(d[i]) - 1
                break
            elif i == len(new_pos) - 1:
                iter += len(atom.GetSymbol())
        counter += 1
    return smi

def canonicalize_repeat_units(repeat_units):
    for i in range(len(repeat_units)):
        for j in range(len(repeat_units[i])):
            if "done" in repeat_units[i][j]:
                repeat_units[i][j] = repeat_units[i][j][0:-4]
                continue
            
            # find and replace descriptors with non-common polymer atom [Bk]
            d = []
            d += re.findall(r"\[.\d+\]", repeat_units[i][j])
            for k in d:
                repeat_units[i][j] = repeat_units[i][j].replace(k,"[Bk]")

            # record the atom indices of where the descriptors are located for canonicalization
            index = []
            m = Chem.MolFromSmiles(repeat_units[i][j])
            for q, atom in enumerate(m.GetAtoms()):
                if atom.GetSymbol() == "Bk":
                    index.append(atom.GetIdx())

            # canonicalize and keep track of the original and final atom indices
            # https://github.com/rdkit/rdkit/issues/794
            # https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg06931.html
            m = Chem.MolFromSmiles(repeat_units[i][j])
            for a in m.GetAtoms():
                a.SetProp("foo", str(a.GetIdx()))
            smi = Chem.MolToSmiles(m)
            order = m.GetPropsAsDict(True,True)["_smilesAtomOutputOrder"]
            canonical_order = list(order)
            repeat_units[i][j] = add_descriptors_back(d, index, canonical_order, smi) 

    return repeat_units

def reorder_repeat_units(repeat_units):
    def replace_nonnested_descriptors(repeat):
        desc_indices = [m.start() for m in re.finditer(r"\[.\d+\]", repeat)]
        nested = get_objects(repeat)
        object_indices = []
        for i in range(len(nested[1])):
            object_indices.append([nested[1][i], nested[1][i] + len(nested[0][i])])
        not_nested = []
        for i in range(len(desc_indices)):
            nested = False
            for j in range(len(object_indices)):
                if desc_indices[i] > object_indices[j][0] and desc_indices[i] < object_indices[j][1]:
                    nested = True
                    break
            if not nested:
                not_nested.append(desc_indices[i])
        for i in range(len(not_nested)):
            repeat = repeat[0:not_nested[i]] + "[Bk]" + repeat[not_nested[i] + 4:]
        repeat = repeat.replace("[Bk]","[]")
        return repeat 
    
    r = []
    for stoch in repeat_units:
        repeat_to_ascii = dict()
        ascii_list = []
        for repeat in stoch:
            ascii = []
            repeat_without_d = replace_nonnested_descriptors(repeat)
            for char in repeat_without_d:
                ascii.append(ord(char))
            ascii_list.append(ascii)
            repeat_to_ascii[repeat] = ascii
            
        # print("ascii list: ", ascii_list)
        ascii_list = sorted(ascii_list)
        ordered_repeats = []
        for repeat in ascii_list:
            for key in repeat_to_ascii:
                if repeat_to_ascii[key] == repeat:
                    ordered_repeats.append(key)
        r.append(ordered_repeats) 

    return r

def relabel_bonding_descriptors(repeat_units, terminals):
    def getCompatible(d):
        if "<" in d:
            return d.replace("<", ">")
        elif ">" in d:
            return d.replace(">", "<")
        else:
            return d        
        
    for i in range(len(repeat_units)):
        original_desc = [terminals[i][0]]
        nested_list = []
        initial_length = []
        for j in range(len(repeat_units[i])):
            initial_length.append(len(repeat_units[i][j]))
            nested = get_objects(repeat_units[i][j])
            nested_list.append(nested)
            # print("before object replacement: ", repeat_units[i][j])
            for n in nested[0]:
                repeat_units[i][j] = repeat_units[i][j].replace(n, "")
            # print("after replacement: ", repeat_units[i][j])
            original_desc += re.findall(r"\[.\d+\]", repeat_units[i][j])

        original_desc += [terminals[i][1]]
        new_desc = copy.deepcopy(original_desc)
        counter = 1
        already_treated = []
        # relabel the bonding id number
        for j in range(1, len(original_desc) - 1):
            if j not in already_treated:
                already_treated.append(j)
                for index in range(len(original_desc)):
                    if original_desc[index] == original_desc[j] or original_desc[index] == getCompatible(original_desc[j]):
                        new_desc[index] = "[" + original_desc[index][1] + str(counter) + "]"
                        already_treated.append(index)
                counter += 1
                
        already_treated = []
        # relabel the bonding descriptors <>$
        for j in range(1, len(new_desc) - 1):
            digit = re.findall(r"\d+", new_desc[j])[0]
            if j not in already_treated:
                already_treated.append(j)
                if ">" in new_desc[j]:
                    already_treated.append(j)
                    new_desc[j] = "[<" + str(digit) + "]"
                    for index in range(len(new_desc)):
                        try:
                            digit2 = re.findall(r"\d+", new_desc[index])[0]
                            if digit2 == digit and index not in already_treated:
                                already_treated.append(index)
                                if "<" in new_desc[index]:
                                    new_desc[index] = "[>" + str(digit2) + "]"
                                else:
                                    new_desc[index] = "[<" + str(digit2) + "]"
                        except:
                            continue
                else:
                    for index in range(len(new_desc)):
                        try:
                            digit2 = re.findall(r"\d+", new_desc[index])[0]
                            if digit2 == digit and index not in already_treated:
                                already_treated.append(index)
                        except:
                            continue
        
        terminals[i][0] = new_desc[0]
        terminals[i][1] = new_desc[-1]
        original_desc = original_desc[1:-1]
        new_desc = new_desc[1:-1]
        counter = 0
        # replace repeat unit descriptors with modified descriptors in new_desc
        for j in range(len(repeat_units[i])):   
            index = 0
            new_repeat = ""
            while index < len(repeat_units[i][j]):
                if index <= len(repeat_units[i][j]) - 4 and counter < len(original_desc) and repeat_units[i][j][index:index + 4] == original_desc[counter]:
                    new_repeat += new_desc[counter]
                    counter += 1
                    index += 4
                else:
                    new_repeat += repeat_units[i][j][index]
                    index += 1
            repeat_units[i][j] = new_repeat

            if len(nested_list[j][0]) > 0:
                index = 0
                q = 0
                counter = 0
                final_string = ""
                while index < len(repeat_units[i][j]):
                    if q < len(nested_list[j][1]) and counter == nested_list[j][1][q]:
                        final_string += nested_list[j][0][q]
                        counter += len(nested_list[j][0][q])
                        q += 1
                    else:
                        final_string += repeat_units[i][j][index]
                        index += 1
                        counter += 1
                repeat_units[i][j] = final_string 
    
    return repeat_units, terminals

def canonicalize_individual_endgroup_config(bigsmiles_or_nested_repeat, repeat_units, terminals):
    bigsmiles_or_nested_repeat = str(bigsmiles_or_nested_repeat)

    ## Replace descriptors with Cf and nested objects with Bk
    bigsmiles_or_nested_repeat = sub_obj_Bk(bigsmiles_or_nested_repeat)[0]
    d = []
    d += re.findall(r"\[.\d+\]", bigsmiles_or_nested_repeat)
    for i in d:
        bigsmiles_or_nested_repeat = bigsmiles_or_nested_repeat.replace(i, "[Cf]")

    ## Record the locations of the descriptors (Cf) and nested objects (Bk) and replace with (Bk)
    object_index = []
    desc_index = []
    m = Chem.MolFromSmiles(bigsmiles_or_nested_repeat)
    for q, atom in enumerate(m.GetAtoms()):
        if atom.GetSymbol() == "Bk":
            object_index.append(atom.GetIdx())
        if atom.GetSymbol() == "Cf":
            desc_index.append(atom.GetIdx())
    bigsmiles_or_nested_repeat = bigsmiles_or_nested_repeat.replace("[Cf]","[Bk]")
    
    # print(bigsmiles_or_nested_repeat)

    # Canonicalize polymer and keep track of how objects (index) are shifted in canonical_order
    m = Chem.MolFromSmiles(bigsmiles_or_nested_repeat)
    for a in m.GetAtoms():
        a.SetProp("foo", str(a.GetIdx()))
    polymer_after = Chem.MolToSmiles(m)
    order = m.GetPropsAsDict(True,True)["_smilesAtomOutputOrder"]
    canonical_order = list(order)

    # print("object: ", repeat_units)
    # print("terminals: ", terminals)
    # print("canonical order: ", canonical_order)
    # print("index: ", index)
    # print("d_index: ", d_index)

    # Determine if during canonicalization, the end groups "flipped"
    forward = []
    for i in object_index:
        if i > 0:
            left_index = canonical_order.index(i - 1)
            if left_index < canonical_order.index(i):
                forward.append(True)
            else:
                forward.append(False)
        elif i < len(canonical_order) - 1:
            right_index = canonical_order.index(i + 1)
            if canonical_order.index(i + 1) > canonical_order.index(i):
                forward.append(True)
            else:
                forward.append(False)
        else:
            forward.append(True)
    
    # add_in relates the pre-canonicalized index to the string
    add_in = dict()
    for i in range(len(object_index)):
        object = "{"
        if forward[i]:
            object += terminals[i][0]
            for k in repeat_units[i]:
                object = object + k + ","
            object = object[0:-1] + terminals[i][1] + "}"
        else:
            object += terminals[i][1]
            for k in repeat_units[i]:
                object = object + k + ","
            object = object[0:-1] + terminals[i][0] + "}"
        add_in[object_index[i]] = object
    for i in range(len(desc_index)):
        add_in[desc_index[i]] = d[i]
    #print("add_in", add_in)

    # canonical_order_dict relates the canonicalized index --> pre-canonicalized index
    canonical_order_dict = dict()
    for i in object_index:
        canonical_order_dict[canonical_order.index(i)] = i
    for i in desc_index:
        canonical_order_dict[canonical_order.index(i)] = i
    canonical_order_dict = collections.OrderedDict(sorted(canonical_order_dict.items()))
    stop = list(canonical_order_dict.keys())
    #print("canonical_order_dict",canonical_order_dict)
    
    # Add objects and descriptors back into canonicalized string
    counter = 0
    iter = 0
    i = 0
    m = Chem.MolFromSmiles(polymer_after)
    for q, atom in enumerate(m.GetAtoms()):
        polymer_after_lower = polymer_after.lower()
        while polymer_after_lower[iter:iter + len(atom.GetSymbol())] != atom.GetSymbol().lower():
            iter += 1
        #print("counter", counter)
        if counter == stop[i]:
            polymer_after = polymer_after[0:iter - 1] + add_in[canonical_order_dict[stop[i]]] + polymer_after[iter + 1 + len(atom.GetSymbol()):]
            iter += len(add_in[canonical_order_dict[stop[i]]]) - 1
            i += 1
            if i == len(stop):
                break
        else:
            iter += len(atom.GetSymbol())
        counter += 1 

    return polymer_after  + "done"

def evaluate_objects(bigsmiles_or_nested_repeat):
    
    def canonicalize(bigsmiles_or_nested_repeat, repeat_units, terminals):
        repeat_units = canonicalize_repeat_units(repeat_units)
        # print(repeat_units)
        repeat_units = reorder_repeat_units(repeat_units)
        # print(repeat_units)
        repeat_units, terminals = relabel_bonding_descriptors(repeat_units, terminals)
        # print(repeat_units, terminals)
        bigsmiles_or_nested_repeat = canonicalize_individual_endgroup_config(bigsmiles_or_nested_repeat, repeat_units, terminals)
        return bigsmiles_or_nested_repeat
        
    repeat_units = []
    terminals = []
    objects = get_objects(bigsmiles_or_nested_repeat)[0]
    for object in objects:
        repeat_in_object = []
        repeats = get_repeats(object)
        for repeat in repeats:
            nested = get_objects(repeat)[0]
            if len(nested) > 0:
                repeat_in_object.append(evaluate_objects(repeat))
            else:
                repeat_in_object.append(repeat)
        repeat_units.append(repeat_in_object)
        left_terminal = object[1:object.find("]")+1]
        right_terminal = object[object.rfind("["):-1]
        terminals.append([left_terminal, right_terminal]) 

    bigsmiles_or_nested_repeat = canonicalize(bigsmiles_or_nested_repeat, repeat_units, terminals)
    return bigsmiles_or_nested_repeat