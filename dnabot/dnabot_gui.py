# -*- coding: utf-8 -*-
"""
Created on Thu May 30 14:35:26 2019

@author: mh2210
"""

import tkinter as tk
from tkinter import filedialog
import sys


class UserDefinedPaths:

    def __init__(self, master, title_name, multiple_files=False):
        master.withdraw()
        if multiple_files:
            master.output = filedialog.askopenfilenames(
                title=title_name, filetypes=(("CSV files", "*.CSV"),
                                             ("all files", "*.*")))
        else:
            master.output = filedialog.askopenfilename(
                title=title_name, filetypes=(("CSV files", "*.CSV"),
                                             ("all files", "*.*")))
        self.output = master.output


class DnabotApp:
    dnabot_title = 'dnabot app'
    intro_text = 'Welcome to the dnabot App! Please follow these instructions:'
    instruction_text1 = '1. Use the dropdown menus below to select wells/columns for ethanol and SOC media, used during purification and transformation, respectively.'
    instruction_text2 = """2. Click the 'GENERATE' button below and on the following window select one csv file detailing constructs. On the subsequent window select up to 6 csv files detailing plates containing parts and linkers.
    If all files are not within one folder, absolute paths should be given."""
    app_font = ("Helvetica", 12)
    trough_wells = ['A{}'.format(x + 1) for x in range(12)]

    def __init__(self, master):
        self.master = master
        self.master.lift()

        # Titles and labels
        self.master.title(DnabotApp.dnabot_title)
        intro = tk.Label(self.master, text=DnabotApp.intro_text,
                         font=DnabotApp.app_font)
        intro.grid(row=0, columnspan=2)
        instruction1 = tk.Label(
            self.master, text=DnabotApp.instruction_text1, font=DnabotApp.app_font)
        instruction1.grid(row=1, columnspan=2, sticky=tk.W)
        instruction2 = tk.Label(
            self.master, text=DnabotApp.instruction_text2, font=DnabotApp.app_font)
        instruction2.grid(row=2, columnspan=2, sticky=tk.W)
        etoh_well_label = tk.Label(
            self.master, text='Trough well for ethanol during purification:', font=DnabotApp.app_font)
        etoh_well_label.grid(row = 3, column = 0, sticky = tk.E)
        soc_column_label=tk.Label(
            self.master, text = 'Deep-well plate column for SOC media during transformation:',
            font = DnabotApp.app_font)
        soc_column_label.grid(row = 4, column = 0, sticky = tk.E)

        # Buttons and menus
        self.quit_status=False
        quit_button=tk.Button(
            self.master, text = 'QUIT', fg = 'red', command = self.quitter, 
            font=DnabotApp.app_font)
        quit_button.grid(row = 5, column = 0)
        generate_button=tk.Button(master, text = 'GENERATE',
                                    command = self.generate, font = DnabotApp.app_font)
        generate_button.grid(row=5, column=1)
        self.etoh_well=tk.StringVar(master)
        self.etoh_well.set(DnabotApp.trough_wells[10])
        etoh_w=tk.OptionMenu(master, self.etoh_well, *tuple(
            DnabotApp.trough_wells[1:11]))
        etoh_w.grid(row=3, column=1, sticky=tk.W)
        etoh_w.config(font=DnabotApp.app_font)
        self.soc_column=tk.StringVar(master)
        self.soc_column.set("1")
        soc_w=tk.OptionMenu(master, self.soc_column, *tuple(
            ['{}'.format(x + 1) for x in range(12)]))
        soc_w.grid(row=4, column=1, sticky=tk.W)
        soc_w.config(font=DnabotApp.app_font)

    def quitter(self):
        self.quit_status=True
        self.master.quit()

    def generate(self):
        self.etoh_well=self.etoh_well.get()
        self.soc_column=self.soc_column.get()
        self.master.quit()


def main():
    #    root = tk.Tk()
    #    construct_path = UserDefinedPaths(root, 'Construct csv file')
    #    root.destroy()
    #    constructs_list = basic_ot2_generator_v11_beta.generate_constructs_list(
    #            construct_path.output)
    #    print(constructs_list)
    root=tk.Tk()
    dnabotinst=DnabotApp(root)
    root.mainloop()
    root.destroy()
    if dnabotinst.quit_status:
        sys.exit("User specified 'QUIT' during app")
    print('Ethanol well is ', dnabotinst.etoh_well)
    print('SOC column is ', dnabotinst.soc_column)


if __name__ == '__main__':
    main()
