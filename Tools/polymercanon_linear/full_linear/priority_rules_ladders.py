from dataclasses import replace
from timeit import repeat
from tkinter import E
import rdkit
from rdkit import Chem
import networkx as nx
import re, copy, collections
import pandas 

import formal_language_non_ladders 
import priority_rules_non_ladders 
from BigSMILES_parser.BigSMILES_BigSmilesObj import BigSMILES

ladder_desc = r"\[(\$|\<|\>)\d+\[(\$|\<|\>)\d+\]\d+\]"

def getCompatible(d):
    if "<" in d:
        return d.replace("<", ">")
    elif ">" in d:
        return d.replace(">", "<")
    else:
        return d

def canonicalize_repeat_units(repeat_units):
    # find and replace descriptors with non-common polymer atom [Bk]
    descriptors = []
    descriptor_locations = [(d.start(0), d.end(0)) for d in re.finditer(ladder_desc, repeat_units)] 
    for d in descriptor_locations:
        descriptors.append(repeat_units[d[0]:d[1]])
    for d in descriptors:
        repeat_units = repeat_units.replace(d,"[Bk]")

    # record the atom indices of where the descriptors are located for canonicalization
    index = []
    m = Chem.MolFromSmiles(repeat_units)
    for q, atom in enumerate(m.GetAtoms()):
        if atom.GetSymbol() == "Bk":
            index.append(atom.GetIdx())

    # canonicalize and keep track of the original and final atom indices
    m = Chem.MolFromSmiles(repeat_units)
    for a in m.GetAtoms():
        a.SetProp("foo", str(a.GetIdx()))
    smi = Chem.MolToSmiles(m)
    order = m.GetPropsAsDict(True,True)["_smilesAtomOutputOrder"]
    canonical_order = list(order)
    repeat_units = priority_rules_non_ladders.add_descriptors_back(descriptors, index, canonical_order, smi)
    
    return repeat_units

def relabel_bonding_descriptors(repeat_unit):
    outer = [r"\[.\d+\[",  r"\[.\d+\]"]
    final = []
    for o in outer:
        d = []
        d += re.findall(o, repeat_unit)
        d2 = copy.deepcopy(d)
        counter = 1
        already_treated = []
        # relabel the bonding id number
        for j in range(len(d)):
            if j not in already_treated:
                already_treated.append(j)
                for k in range(len(d)):
                    if d[k] == d[j] or d[k] == getCompatible(d[j]):
                        d2[k] = "[" + d[k][1] + str(counter) + "]"
                        already_treated.append(k)
                counter += 1
        already_treated = []
        # relabel the bonding descriptors <>$
        for j in range(len(d2)):
            digit = re.findall(r"\d+", d2[j])[0]
            if j not in already_treated:
                already_treated.append(j)
                if "<" in d2[j]:
                    already_treated.append(j)
                    d2[j] = "[>" + str(digit) + "]"
                    for k in range(len(d2)):
                        try:
                            digit2 = re.findall(r"\d+", d2[k])[0]
                            if digit2 == digit and k not in already_treated:
                                already_treated.append(k)
                                if "<" in d2[k]:
                                    d2[k] = "[>" + str(digit2) + "]"
                                else:
                                    d2[k] = "[<" + str(digit2) + "]"
                        except:
                            continue
                else:
                    for k in range(len(d2)):
                        try:
                            digit2 = re.findall(r"\d+", d2[k])[0]
                            if digit2 == digit and k not in already_treated:
                                already_treated.append(k)
                        except:
                            continue
        final.append(d2)
    original = re.findall(r"\[.\d+\[.\d+\]\d+\]", repeat_unit)
    group_id = re.findall(r"\][0-9]\]", repeat_unit)
    f = []
    for i in range(len(final[0])):
        f.append(final[0][i][:-1]+final[1][i]+group_id[i][1:])
    # print(original)
    # print(f)
    index = []
    counter = 0
    for i in range(0, len(repeat_unit)):
        if repeat_unit[i:i+len(f[0])] == original[counter]:
            index.append(i)
            counter += 1
            if counter == 4:
                break
    # print(index)
    # print(repeat_unit)
    repeat_unit = repeat_unit[0:index[0]] + f[0] + \
                    repeat_unit[index[0]+len(f[0]):index[1]] + f[1] + \
                    repeat_unit[index[1]+len(f[1]):index[2]] + f[2] + \
                    repeat_unit[index[2]+len(f[2]):index[3]] + f[3] + \
                    repeat_unit[index[3]+len(f[3]):]
    return repeat_unit

def evaluate_objects(bigsmiles):
    bigsmiles = priority_rules_non_ladders.process(bigsmiles)
    bigsmiles = canonicalize_repeat_units(bigsmiles[1:-1])
    bigsmiles = relabel_bonding_descriptors(bigsmiles)
    return "{" + bigsmiles + "}"
