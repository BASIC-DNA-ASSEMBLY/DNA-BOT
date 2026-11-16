# -*- coding: utf-8 -*-
"""
Created on Thu May 30 17:05:37 2019

@authors: mh2210, gizembuldum, tduigou
"""


def final_well(sample_number: int) -> str:
    """Determines well containing the final sample from sample number."""
    letter = ["A", "B", "C", "D", "E", "F", "G", "H"]
    final_well_column = sample_number // 8 + (1 if sample_number % 8 > 0 else 0)
    final_well_row = letter[sample_number - (final_well_column - 1) * 8 - 1]
    return final_well_row + str(final_well_column)


def final_12wellplate(sample_number):
    """Determines well containing the final sample from sample number for 12 well plate spotting"""
    if sample_number < 1:
        raise ValueError("Sample numbers start at 1.")
    index = (sample_number - 1) % 12
    return "ABC"[index % 3] + str(index // 3 + 1)
