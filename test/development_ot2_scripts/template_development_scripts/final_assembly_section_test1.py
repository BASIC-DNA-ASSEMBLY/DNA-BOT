# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 11:39:43 2019

@author: mh2210
"""

final_assembly_dict = {
    'A1': ['A1', 'B1', 'C1', 'D1'], 'B1': ['E1', 'F1', 'G1', 'H1', 'A2'],
    'C1': ['A2', 'B2', 'C2', 'D2']}

final_assembly_lens = []
for values in final_assembly_dict.values():
    final_assembly_lens.append(len(values))
unique_assemblies_lens = list(set(final_assembly_lens))
master_mix_well_letters = ['A', 'B', 'C', 'D']

x = unique_assemblies_lens[0]
master_mix_well = master_mix_well_letters[(x - 1) // 6] + str(x - 1)
destination_inds = [i for i, lens in enumerate(
    final_assembly_lens) if lens == x]

#print(destination_inds)
#destination_wells = np.array([key for key, value in list(final_assembly_dict.items())])
#print(list(destination_wells[destination_inds]))

for key, values in list(final_assembly_dict.items()):
    print(type(key), values)