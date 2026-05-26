import uproot

file_path = "outputs/treemaker/WbWb/semihad/p8_ee_WW_ecm345.root"

events = uproot.open(file_path)

print(events.keys())

# for key, value in zip(events.keys(), events.values()):
#     print(events[key].keys())


print(events['events;1'].keys())

x = events['events;1']

print(x['jets_R5_p'].array())
#print(x['jet_R5_theta'].array())

