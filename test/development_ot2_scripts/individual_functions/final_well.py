# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 17:15:33 2019

@author: mh2210
"""

def final_well(tray_1_vols):
    """Determines well containing the final sample from tray_1_vols input.

    Args:
    tray_1_vols (list): List of lists containing volumes of each transformation reaction to be spotted. Nested lists correspond with columns.

    """
    sample_num = len([item for sublist in tray_1_vols for item in sublist])
    letter = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    final_well_column = sample_num // 8 + (1 if sample_num % 8 > 0 else 0)
    final_well_row = letter[sample_num - (final_well_column - 1) * 8 - 1]
    return final_well_row + str(final_well_column)

print(final_well([[5, 5, 5, 5, 5, 5, 5] * 12]))
