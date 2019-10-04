# -*- coding: utf-8 -*-
"""
Created on Thu May 30 14:35:26 2019

@author: mh2210
"""

import tkinter as tk
from tkinter import filedialog


class UserInput:

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


root = tk.Tk()
construct_path = UserInput(root, 'Construct csv file')
root.destroy()
root = tk.Tk()
sources_path = UserInput(root, 'Construct csv file', multiple_files=True)
root.destroy()
print(construct_path.output, sources_path.output)
