import uproot
import numpy as np
import matplotlib.pyplot as plt
import awkward as ak
import hist

file_path = "outputs/treemaker/WbWb/semihad/p8_ee_WW_ecm345.root"

events_large = uproot.open(file_path)


# for key, value in zip(events.keys(), events.values()):
#     print(events[key].keys())


print(events_large['events;1'].keys())

events = events_large['events;1']

def tagChecker(events, tag=''):
    print(f'{tag} tag')
    print('Number of jets total ', events['njets_R5'].array())
    print(f'Number of {tag} jets ', events[f'n{tag}jets_R5_true'].array())

    for num in [5, 8, 85, 9]:
        
        print(f'Working point of {num}', events[f'n{tag}jets_R5_WPp{str(num)}'].array())

    #print(x['jet_R5_theta'].array())




# print(events['HardWs_all/HardWs_all.charge'].array())
# print(events['HardWs_all/HardWs_all.generatorStatus'].array())

def check(a):
    # ids = a[:,0]
    # y = ids[ids != 22]
    # length = len(y) > 0
    # print(length)
    # print('HERE')

    b = ak.any(a == 22, axis=-1)
    print(b)


    return()

# check(events['HardWs_all/HardWs_all.generatorStatus'].array())



print(events['WDaughtersHadrons_sel'].array())
print()
print(events['WDaughtersHadrons/WDaughtersHadrons.generatorStatus'].array())

#tagChecker('b')

# print('Number of jets total ', events['njets_R5'].array())
# print(events['jets_R5_pflavor'].array())

# print('D', events['recojet_isD_R5'].array())
# print('U', events['recojet_isU_R5'].array())
# print('S', events['recojet_isS_R5'].array())
# print('C', events['recojet_isC_R5'].array())
# print('B', events['recojet_isB_R5'].array())
# print('G', events['recojet_isG_R5'].array())




# expects input of events and at what level those events are tagged

wp_keys = ['gen', '5', '8', '85', '9']

wp_dict = {
    'gen':'jets_R5_pflavor',
    '5':'bjets_R5_WPp5',
    '8':'bjets_R5_WPp8',
    '85':'bjets_R5_WPp8',
    '9':'bjets_R5_WPp8',
}

x = events['bjets_R5_WPp5'].array()

b_tag_mask = ak.num(events['bjets_R5_WPp5'].array()) > 0


events_btagged = events.arrays(list(wp_dict.values()))[b_tag_mask]




#print(events_btagged['jets_R5_pflavor'])
#print(events_btagged['bjets_R5_WPp5'])

# should i make c++ functions for most of my analysis or is python fine

y = events['recojet_isU_R5']
print(y)






# use gen level info to get events with

def GetCBEventsGen(events, wp):






    f'nbjets_R5_WPp{wp}'


    return()


print(events['recojet_isD_R5'].array())


