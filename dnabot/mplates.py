# -*- coding: utf-8 -*-
"""
Created on Thu May 30 17:05:37 2019

@author: mh2210
"""

def final_well(sample_number):
    """Determines well containing the final sample from sample number.
    
    """
    letter = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    final_well_column = sample_number // 8 + \
        (1 if sample_number % 8 > 0 else 0)
    final_well_row = letter[sample_number - (final_well_column - 1) * 8 - 1]
    return final_well_row + str(final_well_column)