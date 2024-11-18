import pandas as pd
import networkx as nx
from deepchem.feat.smiles_tokenizer import BasicSmilesTokenizer
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') 
import xlwt 
from xlwt import Workbook 

import polymersearch.graphs
import polymersearch.search_tools

def generate_data(to_excel = False):
    database = pd.read_excel("polymerscribe/query-target-dataset.xlsx", "Small Molecule")
    bigsmiles = database["Targets (BigSMILES)"][:400]
    repeats_ab = set()
    repeats_aabb = set()
    repeats_d = set()

    for d in bigsmiles:
        try:
            for o in polymersearch.graphs.get_objects(d)[0]:
                for r in polymersearch.graphs.get_repeats(o):
                    if "<" in r[0] and ">" in r[0]:
                        repeats_ab.add(r[0])
                    elif r[0].count("$") == 2:
                        repeats_d.add(r[0])
                    elif r[0].count("<") == 2 or r[0].count(">") == 2:
                        repeats_aabb.add(r[0])
        except:
            continue
        
    def generate_copolymers():
        copolymers = []
        for i in range(len(repeats_ab)):
            for j in range(i + 1, len(repeats_ab)):
                bigsmiles = "{[]" + repeats_ab[i] + "," + repeats_ab[j] + "[]}"
                copolymers.append(bigsmiles)

        for i in range(len(repeats_ab)):
            for j in range(i + 1, len(repeats_ab)):
                for k in range(j + 1, len(repeats_ab)):
                    bigsmiles = "{[]" + repeats_ab[i] + "," + repeats_ab[j] + "," + repeats_ab[k] + "[]}"
                    copolymers.append(bigsmiles) 

        return copolymers

    repeats_ab = list(repeats_ab)
    copolymers = generate_copolymers()
    print(len(copolymers))

    def generate_blocks():
        diblocks = []
        for i in range(len(repeats_ab)):
            for j in range(i + 1, len(repeats_ab)):
                for k in range(j + 1, len(repeats_ab)):
                    bigsmiles = "{[]" + repeats_ab[i] + "[<]}" + repeats_ab[j][3:-3] + "{[>]" + repeats_ab[k] + "[]}"
                    diblocks.append(bigsmiles) 

        return diblocks

    blocks = generate_blocks()
    print(len(blocks))

    def generate_segmented():
        segmented = []
        for i in range(len(repeats_ab)):
            for j in range(i + 1, len(repeats_ab)):
                for k in range(j + 1, len(repeats_ab)):
                    bigsmiles = "{[]" + repeats_ab[i][:-3] + "{[>]" + repeats_ab[j] + "[<]}" + repeats_ab[k][3:] + "[]}"
                    segmented.append(bigsmiles) 
        return segmented
    
    segmented = generate_segmented()
    print(len(segmented))
    
    def generate_endgroup():
        endgroup = []
        for i in range(len(repeats_ab)):
            for j in range(i + 1, len(repeats_ab)):
                for k in range(j + 1, len(repeats_ab)):
                    bigsmiles = repeats_ab[i][3:-3] + "{[>]" + repeats_ab[j] + "[<]}" + repeats_ab[k][3:-3]
                    endgroup.append(bigsmiles) 
        return endgroup
    
    endgroup = generate_endgroup()
    print(len(endgroup))

    def generate_macrocycle():
        macrocycle = []
        tokenizer = BasicSmilesTokenizer()
        for i in range(len(repeats_ab)):
            for j in range(i + 1, len(repeats_ab)):
                for k in range(j + 1, len(repeats_ab)):
                    tok = tokenizer.tokenize(repeats_ab[i][3:-3])
                    str = ""
                    for index in range(len(tok)):
                        if index == 1:
                            str += "%10"
                        str += tok[index]
                    str = str + "{[>]" + repeats_ab[j] + "[<]}" + repeats_ab[k][3:-3] + "%10"
                    macrocycle.append(str) 
        return macrocycle
    
    macrocycle = generate_macrocycle()
    print(len(macrocycle))   

    def generate_star():
        star = []
        for i in range(len(repeats_ab)):
            repeat = repeats_ab[i][3:-3]
            for j in range(i + 1, len(repeats_ab)):
                for k in range(j + 1, len(repeats_ab)): 
                    arm1 = repeat + "{[>]" + repeats_ab[i] + "[]}"
                    arm2 = repeat + "{[>]" + repeats_ab[j] + "[]}"
                    arm3 = repeat + "{[>]" + repeats_ab[k] + "[]}"
                    core = repeat + "C" + "(" + arm1 + ")" + "(" + arm2 + ")" + "(" + arm3 + ")" 
                    star.append(core)

        return star 
    
    star = generate_star()  
    print(len(star))  

    def generate_graft():
        graft = []
        for i in range(len(repeats_ab)):
            repeat1 = repeats_ab[i][3:-3]
            for j in range(i + 1, len(repeats_ab)):
                for k in range(j + 1, len(repeats_ab)): 
                    arm = repeat1 + "{[>]" + repeats_ab[j] + "[]}"
                    core = "{[]" + repeats_ab[i][:-3] + "C" + "(" + arm + ")" + repeats_ab[k][3:] + "[]}"
                    graft.append(core)
        return graft 

    graft = generate_graft()
    print(len(graft))

    if to_excel:        
        wb = Workbook() 
        
        sheet1 = wb.add_sheet('copolymer') 
        for i in range(len(copolymers)):
            sheet1.write(i, 0, copolymers[i]) 

        sheet2 = wb.add_sheet('blocks') 
        for i in range(len(blocks)):
            sheet2.write(i, 0, blocks[i]) 
        
        sheet3 = wb.add_sheet('segmented') 
        for i in range(len(segmented)):
            sheet3.write(i, 0, segmented[i]) 
        
        sheet4 = wb.add_sheet('endgroup') 
        for i in range(len(endgroup)):
            sheet4.write(i, 0, endgroup[i])

        sheet5 = wb.add_sheet('macrocycle') 
        for i in range(len(macrocycle)):
            sheet5.write(i, 0, macrocycle[i]) 
        
        sheet6 = wb.add_sheet('star') 
        for i in range(len(star)):
            sheet6.write(i, 0, star[i]) 
        
        sheet7 = wb.add_sheet('graft') 
        for i in range(len(graft)):
            sheet7.write(i, 0, graft[i]) 

        wb.save('polymerscribe/polymerscribe_data.xls') 