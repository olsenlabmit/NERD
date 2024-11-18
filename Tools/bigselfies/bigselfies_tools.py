import selfies as sf
import re
import random

from polymersearch.graphs import sub_obj_with_Bk, get_repeats, desc_regex

def insert_bd(repeat, descriptors):
    for d in descriptors:
        repeat = repeat.replace("[Cf]", d, 1)
    return repeat

def insert_object(string, object):
    for o in object:
        string = string.replace("[Bk]", o, 1)
    return string

def bigsmiles_bigselfies(bigsmiles, encoder=True):
    # substitute objects with [Bk]
    # feed to encoder and reinsert stochastic object string
    # iterate through each repeat unit
    # substitute descriptors with [Cf]
    # substitute objects with Bk
    # feed to encoder and reinsert descriptors

    a = sub_obj_with_Bk(bigsmiles)
    bk_sub = a[0]
    stochastic_obj = a[1] 
    
    # SMILES -> SELFIES -> SMILES translation
    try:
        if encoder:
            selfies = sf.encoder(bk_sub) 
        else:
            selfies = sf.decoder(bk_sub) 
        object_list = []
        for repeat_list in stochastic_obj:
            left_terminal = repeat_list[0:repeat_list.find("]") + 1]
            right_terminal = repeat_list[repeat_list.rfind("["):] 
            object = left_terminal
            repeat_list = get_repeats(repeat_list)[0] 
            for repeat in repeat_list:
                descriptor_locations = [(d.start(0), d.end(0)) for d in re.finditer(desc_regex, repeat)] 
                descriptors = [] 
                for d in descriptor_locations:
                    descriptors.append(repeat[d[0]:d[1]])
                for d in descriptors:
                    repeat = repeat.replace(d,"[Cf]") 
                if encoder:
                    print(repeat)
                    repeat_selfies = sf.encoder(repeat) 
                    alphabet=list(sf.split_selfies(repeat_selfies)) # Gets the alphabet of robust symbols
                    print(alphabet) 
                    rnd_selfies=''.join(random.sample(list(alphabet), len(alphabet)))
                    rnd_smiles=sf.decoder(rnd_selfies)
                    print(rnd_smiles)
                else:
                    repeat_selfies = sf.decoder(repeat)
                repeat_selfies = insert_bd(repeat_selfies, descriptors)
                object += repeat_selfies + ","
            object = object[0:-1]
            object += right_terminal
            object_list.append(object)
        bigsmiles = insert_object(selfies, object_list)
        return bigsmiles

    except sf.EncoderError:
        pass  # sf.encoder error! 

    except sf.DecoderError:
        pass  # sf.encoder error! 

