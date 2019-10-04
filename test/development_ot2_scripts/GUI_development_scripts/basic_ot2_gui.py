# -*- coding: utf-8 -*-
"""
Created on Thu May 30 14:35:26 2019

@author: mh2210
"""

import tkinter as tk
from tkinter import filedialog
import basic_ot2_generator_v11_beta

class UserDefinedPaths:

    def __init__(self, master, title_name, multiple_files=False):
        master.withdraw()
        self.title_name = title_name
        self.multiple_files = multiple_files
        if self.multiple_files:
            master.output = filedialog.askopenfilenames(
                title=title_name, filetypes=(("CSV files", "*.CSV"),
                                             ("all files", "*.*")))
        else:
            master.output = filedialog.askopenfilename(
                title=title_name, filetypes=(("CSV files", "*.CSV"),
                                             ("all files", "*.*")))
        self.output = master.output

def main():
    root = tk.Tk()
    construct_path = UserDefinedPaths(root, 'Construct csv file')
    root.destroy()
    constructs_list = basic_ot2_generator_v11_beta.generate_constructs_list(
            construct_path.output)
    print(constructs_list)
    
if __name__ == '__main__':
    main()
