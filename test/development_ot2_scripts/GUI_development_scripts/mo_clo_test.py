# -*- coding: utf-8 -*-
"""
Created on Thu May 30 14:35:26 2019

@author: mh2210
"""

import tkinter as tk
from tkinter import filedialog

def user_input_csv(
        title_name, 
        multiple_files=False):
    """Asks user for input csv filenames.
    
    Args:
        title_name (str): Title given to window
        multiple_files (boolean): Whether multiple files can be chosen by the user
    
    """
    root = tk.Tk()
    if multiple_files:
        root.output = filedialog.askopenfilenames(
                title = title_name, filetypes = (("CSV files","*.CSV"),
                                                 ("all files","*.*")))
    else:
         root.output = filedialog.askopenfilename(
                title = title_name, filetypes = (("CSV files","*.CSV"),
                                                ("all files","*.*")))
    return root.output

construct_path = user_input_csv('Construct csv file')
sources_path = user_input_csv('Sources csv files', multiple_files=True)
print(construct_path, sources_path)
