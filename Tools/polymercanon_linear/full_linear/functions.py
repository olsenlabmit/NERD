from turtle import forward
import pydot
import csv
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
import re
import pandas as pd

from pylex.pylex.automaton import Automaton, AutomatonState
from pylex.pylex.dfa import DFA, DFAState
from pylex.pylex.nfa import NFA, NFAState
import formal_language_non_ladders  
import formal_language_ladders 
import priority_rules_non_ladders 
import priority_rules_ladders 

def visualize_alphabets(alphabet_dictionary, num, flipped):
    visual = []
    legend = []
    k = list(alphabet_dictionary.keys())
    for i in range(len(k)):
        visual.append(Chem.MolFromSmiles(k[i]))
        legend.append(alphabet_dictionary[k[i]])
    img = Draw.MolsToGridImage(visual, molsPerRow=6, subImgSize=(200, 250), useSVG=True, legends=legend)
    filename = "outputs/polymer_#" + str(num) + "_flip_bigsmiles_" + str(int(flipped)) + "_A_alphabets.svg"
    with open(filename, 'w') as f:
        f.write(img.data)

def visualize_automata(g, meaning, num, flipped):
    g = g.print_graphviz()
    text_file = open("outputs/dot.txt", "r")
    data = text_file.read()
    text_file.close()
    data = data[0:data.rindex(";")] + "}"
    graphs = pydot.graph_from_dot_data(data)
    graphs[0].write_svg("outputs/polymer_#" + str(num) + "_flip_bigsmiles_" + str(int(flipped)) + meaning + ".svg")

def priority_rules(bigsmiles):
    contains_ladder_desc = len([(d.start(0), d.end(0)) for d in re.finditer(priority_rules_ladders.ladder_desc, bigsmiles)]) > 0
    
    if contains_ladder_desc: 
        output = priority_rules_non_ladders.process(bigsmiles)
        priority = priority_rules_ladders.evaluate_objects(output)
    else:
        output = priority_rules_non_ladders.process(bigsmiles)
        # if there are implicit endgroups
        if ";" in output:
            priority = priority_rules_non_ladders.implicit_config(output)
        else:
            priority = priority_rules_non_ladders.evaluate_objects(output)[0:-4]
    
    return priority

def transitions_to_automata(initialState, finalState, statesList, alphabet, transition, alphabet_dictionary, i, forward_automata):
    visualize_alphabets(alphabet_dictionary, i, int(forward_automata))

    # generate automata (NFA) from transitions
    states = []
    states_indices = dict()
    for ind in range(len(statesList)):
        states_indices[statesList[ind]] = ind
        states.append(NFAState())
    states[states_indices[finalState]] = NFAState(accepting=True)
    for t_index in transition:
        states[states_indices[t_index[0]]].add_transition(t_index[2], states[states_indices[t_index[1]]])
    g = NFA(states[states_indices[initialState]])

    # minimize NFA --> DFA --> minimal DFA
    visualize_automata(g, "_B_nfa", i, forward_automata)
    g = g.to_dfa()
    visualize_automata(g, "_C_dfa", i, forward_automata)
    g = g.minimized()
    visualize_automata(g, "_D_mdfa", i, forward_automata)

def formal_language(bigsmiles, graph_label = 0):
    contains_ladder_desc = len([(d.start(0), d.end(0)) for d in re.finditer(priority_rules_ladders.ladder_desc, bigsmiles)]) > 0
    
    if contains_ladder_desc > 0:
        initialState, finalState, statesList, alphabet, transition, alphabet_dictionary = formal_language_ladders.generate_transitions(bigsmiles)
        transitions_to_automata(initialState, finalState, statesList, alphabet, transition, alphabet_dictionary, graph_label, forward_automata = 0)
    else:
        # read the BigSMILES from left to right (forward) and right to left (reverse)
        for forward_automata in [False, True]:

            # generate transitions for BigSMILES
            feed = formal_language_non_ladders.process(bigsmiles)
            automata, backbone, o = formal_language_non_ladders.generate_paths(feed, "start", "end", [], dict(), 1)   
            if not forward_automata:
                automata = formal_language_non_ladders.reverse_paths(automata)  
            initialState, finalState, statesList, alphabet, transition, alphabet_dictionary = formal_language_non_ladders.generate_transitions(automata, backbone)
            
            # convert transitions to state machines
            transitions_to_automata(initialState, finalState, statesList, alphabet, transition, alphabet_dictionary, graph_label, forward_automata)

def database_validation():
    paper = [
        [
            "{[][<][Si](CC)(CC)O[>][]}",
            "{[]CC[Si]([<])(CC)O[>][]}",
            "{[][>2][Si](CC)(CC)O[<2][]}"
        ],
        [
            "{[][$]CC=C(C)C[$],[$]CC(C(C)=C)[$],[$]CC(C)(C=C)[$][]}",
            "{[][$]CC(C)=CC[$],[$]C(C(C)=C)C[$],[$]CC(C=C)(C)[$][]}",
            "{[][$]CC(C(C)=C)[$],[$]CC(C)(C=C)[$],[$]CC=C(C)C[$][]}",
            "{[][$]CC=C(C)C[$],[$]CC(C)(C=C)[$],[$]CC(C(C)=C)[$][]}",
            "{[][$2]CC=C(C)C[$2],[$2]CC(C(C)=C)[$2],[$2]CC(C)(C=C)[$2][]}"
        ],
        [
            "{[][<]OCCCCO[<],[>]C(=O)CCCCC(=O)[>][]}",
            "{[]O([<])CCCCO[<],O=C([>])CCCCC([>])=O[]}",
            "{[][>]C(=O)CCCCC(=O)[>],[<]OCCCCO[<][]}",
            "{[][>2]OCCCCO[>2],[<2]C(=O)CCCCC(=O)[<2][]}"
        ],
        [
            "{[][$]CC(C(=O)OCCC)[$][]}",
            "{[]CCCOC(=O)C([$])C[$][]}",
            "{[][$2]CC(C(=O)OCCC)[$2][]}"
        ],
        [
            "{[][$]Cc1c(CC)cc(cc1)C[$][]}",
            "{[][$]Cc5cc(CC)c(cc5)C[$][]}",
            "{[][$5]Cc1c(CC)cc(cc1)C[$5][]}"
        ],
        [
            "{[][$]Cc1ccc(cc1)C[$][]}",
            "{[]C([$])c2ccc(C[$])cc2[]}",
            "{[][$2]Cc1ccc(cc1)C[$2][]}"
        ],
        [
            "{[][<]OCCO[<],[>]C(=O)c1c2ccccc2c(cc1)C(=O)[>][]}",
            "{[][<]OCCO[<],O=C([>])c5c6ccccc6c(cc5)C([>])=O[]}",
            "{[][>]C(=O)c1c2ccccc2c(cc1)C(=O)[>],[<]OCCO[<][]}",
            "{[][>2]OCCO[>2],[<2]C(=O)c1c2ccccc2c(cc1)C(=O)[<2][]}"
        ],
        [
            "{[][$]Cc1c(Cl)cc(cc1)C[$][]}",
            "{[][$]Cc1ccc(C[$])c(Cl)c1[]}",
            "{[][$4]Cc1c(Cl)cc(cc1)C[$4][]}",
        ],
        [
            "{[][$]CC(C)(C(=O)OC1CCCCC1)[$][]}",
            "{[]CC([$])(C(=O)OC5CCCCC5)C[$][]}",
            "{[][$3]CC(C)(C(=O)OC1CCCCC1)[$3][]}"
        ],
        [
            "{[][$]CC(c1c(F)c(F)c(F)c(F)c1F)[$][]}",
            "{[][$]CC([$])c1c(F)c(F)c(F)c(F)c1F[]}",
            "{[][$2]CC(c1c(F)c(F)c(F)c(F)c1F)[$2][]}"
        ],
        [
            "{[][<]Oc1ccc(cc1)Sc2ccc(cc2)O[<],[>]C(=O)[>][]}",
            "{[][<]Oc2ccc(Sc1ccc(O[<])cc1)cc2,[>]C([>])=O[]}",
            "{[][>]C(=O)[>],[<]Oc1ccc(cc1)Sc2ccc(cc2)O[<][]}",
            "{[][<5]Oc1ccc(cc1)Sc2ccc(cc2)O[<5],[>5]C(=O)[>5][]}"
        ],
        [
            "{[][<]Oc1c(C)cc(cc1)C(C)(C)c2cc(Cl)c(cc2)O[<],[>]C(=O)[>][]}",
            "{[]Cc2cc(C(C)(C)c1ccc(O[<])c(Cl)c1)ccc2O[<],C([>])([>])=O[]}",
            "{[][>]C(=O)[>],[<]Oc1c(C)cc(cc1)C(C)(C)c2cc(Cl)c(cc2)O[<][]}",
            "{[][<2]Oc1c(C)cc(cc1)C(C)(C)c2cc(Cl)c(cc2)O[<2],[>2]C(=O)[>2][]}"
        ],
        [
            "{[][<]Oc1ccc(cc1)C2(CCCC2)c3ccc(cc3)O[<],[>]C(=O)[>][]}",
            "{[][<]Oc3ccc(C2(c1ccc(O[<])cc1)CCCC2)cc3,[>]C(=O)[>][]}",
            "{[][>]C(=O)[>],[<]Oc1ccc(cc1)C2(CCCC2)c3ccc(cc3)O[<][]}",
            "{[][<2]Oc1ccc(cc1)C2(CCCC2)c3ccc(cc3)O[<2],[>2]C(=O)[>2][]}",
        ],
        [
            "{[][<]N1C(=O)c2ccc(cc2C1=O)Oc3ccc(cc3)Sc4ccc(cc4)Oc5ccc6C(=O)N(C(=O)c6c5)[<],[>]c1cc(ccc1)C(=O)c2cc(ccc2)[>][]}",
            "{[]O=c6c5ccc(Oc4ccc(Sc3ccc(Oc2ccc1c(=O)n([<])c(=O)c1c2)cc3)cc4)cc5c(=O)n6[<],[>]c(ccc1)cc1C(=O)c(ccc3)cc3[>][]}",
            "{[][>]c1cc(ccc1)C(=O)c2cc(ccc2)[>],[<]N1C(=O)c2ccc(cc2C1=O)Oc3ccc(cc3)Sc4ccc(cc4)Oc5ccc6C(=O)N(C(=O)c6c5)[<][]}",
            "{[][<2]N1C(=O)c2ccc(cc2C1=O)Oc3ccc(cc3)Sc4ccc(cc4)Oc5ccc6C(=O)N(C(=O)c6c5)[<2],[>2]c1cc(ccc1)C(=O)c2cc(ccc2)[>2][]}",
        ],
        [
            "{[][<]N1C(=O)c2cc(ccc2C1=O)C(C(F)(F)F)(C(F)(F)F)c3cc4C(=O)N(C(=O)c4cc3)[<],[>]c1cc(ccc1)C(=O)c2cc(ccc2)[>][]}",
            "{[]O=c4c3ccc(C(c2ccc1c(=O)n([<])c(=O)c1c2)(C(F)(F)F)C(F)(F)F)cc3c(=O)n4[<],[>]c2cccc(C(=O)c1cccc([>])c1)c2[]}",
            "{[][>]c1cc(ccc1)C(=O)c2cc(ccc2)[>],[<]N1C(=O)c2cc(ccc2C1=O)C(C(F)(F)F)(C(F)(F)F)c3cc4C(=O)N(C(=O)c4cc3)[<][]}",
            "{[][<3]N1C(=O)c2cc(ccc2C1=O)C(C(F)(F)F)(C(F)(F)F)c3cc4C(=O)N(C(=O)c4cc3)[<3],[>3]c1cc(ccc1)C(=O)c2cc(ccc2)[>3][]}",
        ],
        [
            "{[][<]c1cc(c2ccccc2)c3cc(ccc3n1)Oc4cc5c(c6ccccc6)cc(nc5cc4)[<],[>]c1ccc(cc1)[>][]}",
            "{[][<]c6cc(c1ccccc1)c5cc(Oc4ccc3nc([<])cc(c2ccccc2)c3c4)ccc5n6,[>]c1ccc([>])cc1[]}",
            "{[][>]c1ccc(cc1)[>],[<]c1cc(c2ccccc2)c3cc(ccc3n1)Oc4cc5c(c6ccccc6)cc(nc5cc4)[<][]}",
            "{[][<4]c1cc(c2ccccc2)c3cc(ccc3n1)Oc4cc5c(c6ccccc6)cc(nc5cc4)[<4],[>4]c1ccc(cc1)[>4][]}",
        ],
        [
            "{[][<]C(=O)c1ccc(cc1)C2(OC(=O)c3ccccc23)c4ccc(cc4)C(=O)[<],[>]Oc1ccc(cc1)C2(OC(=O)c3ccccc23)c4ccc(cc4)O[>][]}",
            "{[][<]C(=O)c4ccc(C3(c1ccc(C([<])=O)cc1)OC(=O)c2ccccc23)cc4,[>]Oc4ccc(C3(c1ccc(O[>])cc1)OC(=O)c2ccccc23)cc4[]}",
            "{[][>]Oc1ccc(cc1)C2(OC(=O)c3ccccc23)c4ccc(cc4)O[>],[<]C(=O)c1ccc(cc1)C2(OC(=O)c3ccccc23)c4ccc(cc4)C(=O)[<][]}",
            "{[][>2]C(=O)c1ccc(cc1)C2(OC(=O)c3ccccc23)c4ccc(cc4)C(=O)[>2],[<2]Oc1ccc(cc1)C2(OC(=O)c3ccccc23)c4ccc(cc4)O[<2][]}"
        ],
        [
            "{[][<]c1c(c2ccccc2)c(c3ccccc3)c4cc(ccc4n1)Oc5cc6c(c7ccccc7)c(c8ccccc8)c(nc6cc5)[<],[>]c1ccc(cc1)[>][]}",
            "{[][<]c7nc5ccc(Oc4ccc3nc([<])c(c1ccccc1)c(c2ccccc2)c3c4)cc5c(c6ccccc6)c7c8ccccc8,[>]c1ccc([>])cc1[]}",
            "{[][>]c1ccc(cc1)[>],[<]c1c(c2ccccc2)c(c3ccccc3)c4cc(ccc4n1)Oc5cc6c(c7ccccc7)c(c8ccccc8)c(nc6cc5)[<][]}",
            "{[][>3]c1c(c2ccccc2)c(c3ccccc3)c4cc(ccc4n1)Oc5cc6c(c7ccccc7)c(c8ccccc8)c(nc6cc5)[>3],[<3]c1ccc(cc1)[<3][]}",
        ],
        [
            "{[][<]c1c(c2ccccc2)c(c3ccccc3)c4cc(ccc4n1)Oc5cc6c(c7ccccc7)c(c8ccccc8)c(nc6cc5)[<],[>]c1ccc(cc1)c2ccc(cc2)[>][]}",
            "{[][<]c7nc5ccc(Oc4ccc3nc([<])c(c1ccccc1)c(c2ccccc2)c3c4)cc5c(c6ccccc6)c7c8ccccc8,[>]c2ccc(c1ccc([>])cc1)cc2[]}",
            "{[][>]c1ccc(cc1)c2ccc(cc2)[>],[<]c1c(c2ccccc2)c(c3ccccc3)c4cc(ccc4n1)Oc5cc6c(c7ccccc7)c(c8ccccc8)c(nc6cc5)[<][]}",
            "{[][>8]c1c(c2ccccc2)c(c3ccccc3)c4cc(ccc4n1)Oc5cc6c(c7ccccc7)c(c8ccccc8)c(nc6cc5)[>8],[<8]c1ccc(cc1)c2ccc(cc2)[<8][]}"
        ],
        [
            "{[][<]c1nc2cc(ccc2nc1)c3cc4nc(cnc4cc3)[<],[>]c1ccc(cc1)[>][]}",
            "{[][<]c4cnc3ccc(c2ccc1ncc([<])nc1c2)cc3n4,[>]c1ccc([>])cc1[]}",
            "{[][<]c1nc2cc(ccc2nc1)c3cc4nc(cnc4cc3)[<],[>]c1ccc(cc1)[>][]}",
            "{[][>2]c1ccc(cc1)[>2],[<2]c1nc2cc(ccc2nc1)c3cc4nc(cnc4cc3)[<2][]}"
        ],
        [
            "[H]{[>][<]OC(=O)OC(C)C[>][<]}OCC{[$][$]CC(c1ccccc1)[$][$]}C(C)CC",
            "[H]{[>]CC(C[>])OC(=O)O[<][<]}OCC{[$][$]CC([$])c1ccccc1[$]}C(C)CC",
            "[H]{[>3][<3]OC(=O)OC(C)C[>3][<3]}OCC{[$2][$2]CC(c1ccccc1)[$2][$2]}C(C)CC",
            "CCC(C){[$][$]CC(c1ccccc1)[$][$]}CCO{[>][<]CC(C)OC(=O)O[>][<]}[H]"
        ],
        [
            "CCC(C){[$][$]CC(C1CCCCC1)[$][$]}{[$][$]CCCC[$],[$]CC(CC)[$][$]}[H]",
            "CCC(C){[$][$]CC([$])C1CCCCC1[$]}{[$][$]CCCC[$],CCC([$])C[$][$]}[H]",
            "CCC(C){[$][$]CC(C1CCCCC1)[$][$]}{[$][$]CC(CC)[$],[$]CCCC[$][$]}[H]",
            "CCC(C){[$2][$2]CC(C1CCCCC1)[$2][$2]}{[$3][$3]CCCC[$3],[$3]CC(CC)[$3][$3]}[H]",
            "[H]{[$][$]CCCC[$],[$]CC(CC)[$][$]}{[$][$]CC(C1CCCCC1)[$][$]}C(C)CC"
        ],
        [
            "CCC(C){[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}{[>][<][Si](C)(C)O[>][<]}[Si](C)(C)C",
            "CCC(C){[$][$]CC(C)CC[$],CC(C)C([$])C[$],CCC(C)([$])C[$][$]}{[>]C[Si](C)([<])O[>][<]}[Si](C)(C)C",
            "CCC(C){[$][$]CCC(C)C[$],[$]CC(C)(CC)[$],[$]CC(C(C)C)[$][$]}{[>][<][Si](C)(C)O[>][<]}[Si](C)(C)C",
            "CCC(C){[$3][$3]CCC(C)C[$3],[$3]CC(C(C)C)[$3],[$3]CC(C)(CC)[$3][$3]}{[>2][<2][Si](C)(C)O[>2][<2]}[Si](C)(C)C",
            "C[Si](C)(C){[<][<][Si](C)(C)O[>][>]}{[$][$]CCC(C)C[$],[$]CC(C(C)C)[$],[$]CC(C)(CC)[$][$]}C(C)CC"
        ],
        [
            "[H]{[>][<]OC(C)C(=O)OC(C)C(=O)[>][<]}OCC{[$][$]C\C=C(C)/C[$],[$]C\C=C(C)\C[$],[$]CC(C(C)=C)[$],[$]CC(C)(C=C)[$][$]}C(C)CC",
            "[H]{[>]CC(O[<])C(=O)OC(C)C([>])=O[<]}OCC{[$][$]C\C=C(C)/C[$],[$]C\C=C(C)\C[$],[$]CC(C(C)=C)[$],[$]CC(C)(C=C)[$][$]}C(C)CC",
            "[H]{[>][<]OC(C)C(=O)OC(C)C(=O)[>][<]}OCC{[$][$]C\C=C(C)\C[$],[$]C\C=C(C)/C[$],[$]CC(C)(C=C)[$],[$]CC(C(C)=C)[$][$]}C(C)CC",
            "[H]{[>4][<4]OC(C)C(=O)OC(C)C(=O)[>4][<4]}OCC{[$2][$2]C\C=C(C)/C[$2],[$2]C\C=C(C)\C[$2],[$2]CC(C(C)=C)[$2],[$2]CC(C)(C=C)[$2][$2]}C(C)CC",
            "CCC(C){[$][$]C\C=C(C)/C[$],[$]C\C=C(C)\C[$],[$]CC(C(C)=C)[$],[$]CC(C)(C=C)[$][$]}CCO{[>][>]OC(C)C(=O)OC(C)C(=O)[<][<]}[H]"
        ],
        [
            "{[][$]CC(C)(C(=O)OC)[$],[$]CC(C)(C(=O)OCCOC(=O)C(C)(C){[$][$]CC(c1ccccc1)[$][]})[$][]}",
            "{[]COC(=O)C(C)([$])C[$],[$]CC(C)(C(=O)OCCOC(=O)C(C)(C){[$][$]C(c5ccccc5)C[$][]})[$][]}",
            "{[][$]CC(C)(C(=O)OCCOC(=O)C(C)(C){[$][$]CC(c1ccccc1)[$][]})[$],[$]CC(C)(C(=O)OC)[$][]}",
            "{[][$3]CC(C)(C(=O)OC)[$3],[$3]CC(C)(C(=O)OCCOC(=O)C(C)(C){[$2][$2]CC(c1ccccc1)[$2][]})[$3][]}",
            "{[][$]CC(C)(C(=O)OC)[$],{[][$]CC(c1ccccc1)[$][$]}C(C)(C)C(=O)OCCOC(=O)C(C)([$])C[$][]}"
        ],
        [
            "COCCO{[>][<]C(O{[<][>]CCO[<][>]}C)CCCCCO[>],[<]C(=O)CCCCCO[>][<]}",
            "COCCO{[>]C{[<][<]CCO[>][>]}OC([<])CCCCCO[>],[>]OCCCCCC([<])=O[<]}",
            "COCCO{[>][<]C(=O)CCCCCO[>],[<]C(O{[<][>]CCO[<][>]}C)CCCCCO[>][<]}",
            "COCCO{[>2][<2]C(O{[<][>]CCO[<][>]}C)CCCCCO[>2],[<2]C(=O)CCCCCO[>2][<2]}",
            "{[<]C{[<][<]CCO[>][>]}OC([<])CCCCCO[>],[<]C(=O)CCCCCO[>][>]}OCCOC"
        ],
        [
            "{[][$]CC(c(cc1)ccc1CCCC{[>][<]OCC[>][<]}N)[$],[$]CC(c1ccccc1)[$][]}",
            "{[][$]CC(c(cc1)ccc1CCCC{[>][>]CCO[<][<]}N)[$],[$]CC([$])c1ccccc1[]}",
            "{[][$]CC(c1ccccc1)[$],[$]CC(c(cc1)ccc1CCCC{[>][<]OCC[>][<]}N)[$][]}",
            "{[][$2]CC(c(cc1)ccc1CCCC{[>3][<3]OCC[>3][<3]}N)[$2],[$2]CC(c1ccccc1)[$2][]}",
            "{[]N{[>][<]CCO[>][<]}CCCCc1ccc(C([$])C[$])cc1,[$]CC(c1ccccc1)[$][]}"
        ],
        [
            "{[][$]CC(c(cc1)ccc1CCCC{[>][<]OCC[>][<]}C(=O)O)[$],[$]CC(c1ccccc1)[$][]}",
            "{[][$]CC(c(cc5)ccc5CCCC{[>][>]CCO[<][<]}C(O)=O)[$],c5ccccc5C([$])C[$][]}",
            "{[][$]CC(c1ccccc1)[$],[$]CC(c(cc1)ccc1CCCC{[>][<]OCC[>][<]}C(=O)O)[$][]}",
            "{[][$5]CC(c(cc1)ccc1CCCC{[>][<]OCC[>][<]}C(=O)O)[$5],[$5]CC(c1ccccc1)[$5][]}",
            "{[]O=C(O){[>][<]CCO[>][<]}CCCCc1ccc(C([$])C[$])cc1,[$]CC(c1ccccc1)[$][]}",
            "{[]O=C(O){[>][>]OCC[<][<]}CCCCc1ccc(C([$])C[$])cc1,[$]CC(c1ccccc1)[$][]}"
        ],
        [
            "{[][<]{[>][<]C(=O)Nc(cc1)ccc1Cc(cc2)ccc2NC(=O)[<],[>]O{[<][>]CC(C)O[<][>]}[>][>]}[<],[>]NCCN[>][]}",
            "{[][<]{[>][<]C(=O)Nc(cc5)ccc5Cc(cc6)ccc6NC(=O)[<],[>]{[<][>]OC(C)C[<][>]}O[>][>]}[<],[>]NCCN[>][]}",
            "{[][>]NCCN[>],[<]{[>][>]O{[<][>]CC(C)O[<][>]}[>],[<]C(=O)Nc(cc1)ccc1Cc(cc2)ccc2NC(=O)[<][>]}[<][]}",
            "{[][<3]{[>2][<2]C(=O)Nc(cc1)ccc1Cc(cc2)ccc2NC(=O)[<2],[>2]O{[<5][>5]CC(C)O[<5][>5]}[>2][>2]}[<3],[>3]NCCN[>3][]}",
            "{[][<]{[>][<]C(=O)Nc(cc1)ccc1Cc(cc2)ccc2NC(=O)[<],[>]O{[<][<]OC(C)C[>][>]}[>][>]}[<],[>]NCCN[>][]}"
        ],
        [
            "{[][<]C(=O)Nc(cc1)ccc1Cc(cc2)ccc2NC(=O)[<],[>]NCCC{[>][<][Si](C)(C)O[>][<]}[Si](C)(C)CCCN[>][]}",
            "{[]O=C([<])Nc5ccc(Cc1ccc(NC(=O)[<])cc1)cc5,[>]NCCC[Si](C)(C){[>][<]O[Si](C)(C)[>][<]}CCCN[>][]}",
            "{[][>]NCCC{[>][<][Si](C)(C)O[>][<]}[Si](C)(C)CCCN[>],[<]C(=O)Nc(cc1)ccc1Cc(cc2)ccc2NC(=O)[<][]}",
            "{[][<2]C(=O)Nc(cc1)ccc1Cc(cc2)ccc2NC(=O)[<2],[>2]NCCC{[>3][<3][Si](C)(C)O[>3][<3]}[Si](C)(C)CCCN[>2][]}",
            "{[][<]C(=O)Nc(cc1)ccc1Cc(cc2)ccc2NC(=O)[<],[>]NCCC{[>][>]O[Si](C)(C)[<][<]}[Si](C)(C)CCCN[>][]}"
        ],
        [
            "{[][<]NNC(=O)C(=O)NCCC[Si](C)(C){[<][>]O[Si](C)(C)[<][>]}CCCNC(=O)C(=O)NN[<],[>]C(=O)NC(CC1)CCC1CC(CC2)CCC2NC(=O)[>][]}",
            "{[][<]NNC(=O)C(=O)NCCC{[<][>][Si](C)(C)O[<][>]}[Si](C)(C)CCCNC(=O)C(=O)NN[<],[>]C(=O)NC(CC3)CCC3CC(CC4)CCC4NC(=O)[>][]}",
            "{[][>]C(=O)NC(CC1)CCC1CC(CC2)CCC2NC(=O)[>],[<]NNC(=O)C(=O)NCCC[Si](C)(C){[<][>]O[Si](C)(C)[<][>]}CCCNC(=O)C(=O)NN[<][]}",
            "{[][<3]NNC(=O)C(=O)NCCC[Si](C)(C){[<2][>2]O[Si](C)(C)[<2][>2]}CCCNC(=O)C(=O)NN[<3],[>3]C(=O)NC(CC1)CCC1CC(CC2)CCC2NC(=O)[>3][]}",
            "{[][<]NNC(=O)C(=O)NCCC[Si](C)(C){[<][<][Si](C)(C)O[>][>]}CCCNC(=O)C(=O)NN[<],[>]C(=O)NC(CC1)CCC1CC(CC2)CCC2NC(=O)[>][]}"
        ],
        [
            "{[][<]NCCC{[>][<][Si](C)(C)O[>][<]}[Si](C)(C)CCCN[<],[>]C(=O)NC(CC1)CCC1CC(CC2)CCC2NC(=O)[>][]}",
            "{[][<]NCCC[Si](C)(C){[>][<]O[Si](C)(C)[>][<]}CCCN[<],O=C([>])NC2CCC(CC1CCC(NC(=O)[>])CC1)CC2[]}",
            "{[][>]C(=O)NC(CC1)CCC1CC(CC2)CCC2NC(=O)[>],[<]NCCC{[>][<][Si](C)(C)O[>][<]}[Si](C)(C)CCCN[<][]}",
            "{[][<4]NCCC{[>2][<2][Si](C)(C)O[>2][<2]}[Si](C)(C)CCCN[<4],[>4]C(=O)NC(CC1)CCC1CC(CC2)CCC2NC(=O)[>4][]}",
            "{[][<]NCCC{[>][>]O[Si](C)(C)[<][<]}[Si](C)(C)CCCN[<],[>]C(=O)NC(CC1)CCC1CC(CC2)CCC2NC(=O)[>][]}"
        ],
        [
            "{Cc1c([<1[<1]1])c([<1[<1]1])c(C)c3c1C4C2CC(C([>1[>1]2])C2[>1[>1]2])C34}",
            "{Cc1c([<1[<1]1])c([<1[<1]1])c(C)c2c1C1C3CC(C([>1[>1]2])C3[>1[>1]2])C21}", #SMILES
            "{Cc1c([<1[>2]1])c([<1[>2]1])c(C)c3c1C4C2CC(C([>1[<2]2])C2[>1[<2]2])C34}", # inner relabeled
            "{Cc1c([<5[<1]1])c([<5[<1]1])c(C)c3c1C4C2CC(C([>5[>1]2])C2[>5[>1]2])C34}", # outer relabeled
        ],
        [
            "{Cc1c([<1[<1]1])c([<1[<1]1])c(N(C)C)c3c1C4C2CC(C([>1[>1]2])C2[>1[>1]2])C34}",
            "{Cc1c([<1[<1]1])c([<1[<1]1])c(N(C)C)c2c1C1C3CC(C([>1[>1]2])C3[>1[>1]2])C21}", #SMILES
            "{Cc1c([<1[<6]1])c([<1[<6]1])c(N(C)C)c3c1C4C2CC(C([>1[>6]2])C2[>1[>6]2])C34}", #inner relabeled 
            "{Cc1c([>1[<1]1])c([>1[<1]1])c(N(C)C)c3c1C4C2CC(C([<1[>1]2])C2[<1[>1]2])C34}" # outer relabeled
        ],
        [
            "{CCCCCCCCOc1c([<1[<1]1])c([<1[<1]1])c(OCCCCCCCC)c3c1C4C2CC(C([>1[>1]2])C2[>1[>1]2])C34}",
            "{CCCCCCCCOc1c([<1[<1]1])c([<1[<1]1])c(OCCCCCCCC)c2c1C1C3CC(C([>1[>1]2])C3[>1[>1]2])C21}", # SMILES
            "{CCCCCCCCOc1c([<1[<3]1])c([<1[<3]1])c(OCCCCCCCC)c3c1C4C2CC(C([>1[>3]2])C2[>1[>3]2])C34}", # inner relabeled
            "{CCCCCCCCOc1c([<7[<1]1])c([<7[<1]1])c(OCCCCCCCC)c3c1C4C2CC(C([>7[>1]2])C2[>7[>1]2])C34}" # outer relabeled
        ],
        [
            "{Cc1cc2c(c([<1[<1]1])c1[<1[<1]1])C3CC2C([>1[>1]2])C3[>1[>1]2]}",
            "{Cc1cc2c(c([<1[<1]1])c1[<1[<1]1])C1CC2C([>1[>1]2])C1[>1[>1]2]}", # SMILES
            "{Cc1cc2c(c([<1[>4]1])c1[<1[>4]1])C3CC2C([>1[<4]2])C3[>1[<4]2]}", # inner relabeled
            "{Cc1cc2c(c([<2[<1]1])c1[<2[<1]1])C3CC2C([>2[>1]2])C3[>2[>1]2]}" # outer relabeled
        ]
    ]

    count = 0
    for polymer_inputs in paper:
        count += 1
        output = []
        for p in polymer_inputs:
            output.append(priority_rules(p))
        print(count, output[0], len(set(output)) == 1)
        # print(set(output))

if __name__ == "__main__":
    # function runs through 36 examples through priority rules canonicalization
    database_validation()
    quit() 
    bigsmiles = ["{Cc1c([>1[>1]1])c([>1[>1]1])c(N(C)C)c2c1C1C3CC(C([<1[<1]2])C3[<1[<1]2])C21}"]
    for i in range(len(bigsmiles)):
        formal_language(bigsmiles[i], i)
        # priority = priority_rules(bigsmiles[i])
        # print(priority)