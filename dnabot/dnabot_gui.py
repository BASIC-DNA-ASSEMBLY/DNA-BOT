# -*- coding: utf-8 -*-
"""
Created on Thu May 30 14:35:26 2019

@author: mh2210
"""

import tkinter as tk
from tkinter import filedialog, ttk
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
    dnabot_title = "DNA-BOT Protocol Generator"
    intro_text = "Welcome to the DNA-BOT Protocol Generator!"
    instruction_text1 = '1. Select the ethanol well for purification (2_purification script).'
    instruction_text2 = """2. Enter the paths to your constructs and source CSV files containing BASIC parts and linkers."""
    
    # Modern color scheme
    colors = {
        'primary': '#2c3e50',      # Dark blue-gray
        'secondary': '#3498db',    # Blue
        'accent': '#e74c3c',       # Red
        'success': '#27ae60',      # Green
        'background': '#ecf0f1',   # Light gray
        'text': '#2c3e50',         # Dark text
        'light_text': '#7f8c8d'    # Light text
    }
    
    app_font = ("Segoe UI", 10)
    title_font = ("Segoe UI", 16, "bold")
    header_font = ("Segoe UI", 12, "bold")
    trough_wells = ['A{}'.format(x + 1) for x in range(12)]

    def __init__(self, master):
        self.master = master
        self.master.lift()
        
        # Configure window
        self.master.title(DnabotApp.dnabot_title)
        self.master.configure(bg=DnabotApp.colors['background'])
        self.master.geometry("600x500")
        self.master.resizable(True, True)
        
        # Configure grid weights for responsive layout
        self.master.grid_columnconfigure(0, weight=1)
        self.master.grid_columnconfigure(1, weight=2)
        
        # Create main frame with padding
        main_frame = tk.Frame(self.master, bg=DnabotApp.colors['background'])
        main_frame.grid(row=0, column=0, columnspan=2, sticky="nsew", padx=20, pady=20)
        main_frame.grid_columnconfigure(0, weight=1)
        main_frame.grid_columnconfigure(1, weight=2)

        # Title banner
        title_frame = tk.Frame(main_frame, bg=DnabotApp.colors['primary'], height=60)
        title_frame.grid(row=0, column=0, columnspan=2, sticky="ew", pady=(0, 20))
        title_frame.grid_propagate(False)
        
        title_label = tk.Label(title_frame, text=DnabotApp.dnabot_title, 
                              font=DnabotApp.title_font, 
                              fg='white', bg=DnabotApp.colors['primary'])
        title_label.pack(expand=True)
        
        # Welcome text
        intro = tk.Label(main_frame, text=DnabotApp.intro_text,
                         font=DnabotApp.header_font, 
                         fg=DnabotApp.colors['text'],
                         bg=DnabotApp.colors['background'],
                         wraplength=550)
        intro.grid(row=1, column=0, columnspan=2, pady=(0, 15), sticky="ew")
        
        # Instructions section
        instructions_frame = tk.Frame(main_frame, bg=DnabotApp.colors['background'])
        instructions_frame.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(0, 20))
        
        instruction1 = tk.Label(instructions_frame, text=DnabotApp.instruction_text1, 
                               font=DnabotApp.app_font,
                               fg=DnabotApp.colors['text'],
                               bg=DnabotApp.colors['background'],
                               wraplength=550)
        instruction1.pack(anchor="w")
        
        instruction2 = tk.Label(instructions_frame, text=DnabotApp.instruction_text2, 
                               font=DnabotApp.app_font,
                               fg=DnabotApp.colors['text'],
                               bg=DnabotApp.colors['background'],
                               wraplength=550)
        instruction2.pack(anchor="w", pady=(5, 0))
        
        # Settings section
        settings_frame = tk.Frame(main_frame, bg=DnabotApp.colors['background'])
        settings_frame.grid(row=3, column=0, columnspan=2, sticky="ew", pady=(0, 20))
        
        # Ethanol well selection
        etoh_well_label = tk.Label(settings_frame, text='Ethanol well for purification:', 
                                  font=DnabotApp.app_font,
                                  fg=DnabotApp.colors['text'],
                                  bg=DnabotApp.colors['background'])
        etoh_well_label.grid(row=0, column=0, sticky=tk.E, padx=(0, 10), pady=5)
        
        self.etoh_well = tk.StringVar()
        self.etoh_well.set(DnabotApp.trough_wells[10])
        etoh_w = ttk.Combobox(settings_frame, textvariable=self.etoh_well, 
                              values=DnabotApp.trough_wells[1:11],
                              font=DnabotApp.app_font, state="readonly", width=15)
        etoh_w.grid(row=0, column=1, sticky=tk.W, pady=5)
        
        # Water well selection
        water_well_label = tk.Label(settings_frame, text='Water well for purification:', 
                                   font=DnabotApp.app_font,
                                   fg=DnabotApp.colors['text'],
                                   bg=DnabotApp.colors['background'])
        water_well_label.grid(row=0, column=2, sticky=tk.E, padx=(20, 10), pady=5)
        
        self.water_well = tk.StringVar()
        self.water_well.set(DnabotApp.trough_wells[11])
        water_w = ttk.Combobox(settings_frame, textvariable=self.water_well, 
                               values=DnabotApp.trough_wells[1:11],
                               font=DnabotApp.app_font, state="readonly", width=15)
        water_w.grid(row=0, column=3, sticky=tk.W, pady=5)
        
        # Thermocycler generation selection
        thermocycler_label = tk.Label(settings_frame, text='Thermocycler generation:', 
                                     font=DnabotApp.app_font,
                                     fg=DnabotApp.colors['text'],
                                     bg=DnabotApp.colors['background'])
        thermocycler_label.grid(row=1, column=0, sticky=tk.E, padx=(0, 10), pady=5)
        
        self.thermocycler_gen = tk.StringVar()
        self.thermocycler_gen.set("GEN2")
        tc_w = ttk.Combobox(settings_frame, textvariable=self.thermocycler_gen, 
                            values=['None', 'GEN1', 'GEN2'],
                            font=DnabotApp.app_font, state="readonly", width=15)
        tc_w.grid(row=1, column=1, sticky=tk.W, pady=5)
        
        # File inputs section
        files_frame = tk.Frame(main_frame, bg=DnabotApp.colors['background'])
        files_frame.grid(row=4, column=0, columnspan=2, sticky="ew", pady=(0, 20))
        
        # Constructs file input
        construct_label = tk.Label(files_frame, text='Constructs CSV file path:', 
                                 font=DnabotApp.app_font,
                                 fg=DnabotApp.colors['text'],
                                 bg=DnabotApp.colors['background'])
        construct_label.grid(row=0, column=0, sticky=tk.E, padx=(0, 10), pady=5)
        
        self.construct_path = tk.StringVar()
        construct_entry = tk.Entry(files_frame, textvariable=self.construct_path, 
                                 width=40, font=DnabotApp.app_font,
                                 relief="solid", bd=1)
        construct_entry.grid(row=0, column=1, sticky="ew", pady=5)
        
        construct_browse = tk.Button(files_frame, text='Browse', 
                                   fg='white', bg=DnabotApp.colors['secondary'],
                                   command=self.browse_constructs, 
                                   font=DnabotApp.app_font,
                                   relief="flat", bd=0,
                                   padx=15, pady=2,
                                   activebackground='#2980b9',
                                   activeforeground='white')
        construct_browse.grid(row=0, column=2, sticky="w", padx=(10, 0), pady=5)
        
        # Sources files input
        sources_label = tk.Label(files_frame, text='Source CSV files (comma-separated):', 
                               font=DnabotApp.app_font,
                               fg=DnabotApp.colors['text'],
                               bg=DnabotApp.colors['background'])
        sources_label.grid(row=1, column=0, sticky=tk.E, padx=(0, 10), pady=5)
        
        self.sources_paths = tk.StringVar()
        sources_entry = tk.Entry(files_frame, textvariable=self.sources_paths, 
                               width=40, font=DnabotApp.app_font,
                               relief="solid", bd=1)
        sources_entry.grid(row=1, column=1, sticky="ew", pady=5)
        
        sources_browse = tk.Button(files_frame, text='Browse', 
                                 fg='white', bg=DnabotApp.colors['secondary'],
                                 command=self.browse_sources, 
                                 font=DnabotApp.app_font,
                                 relief="flat", bd=0,
                                 padx=15, pady=2,
                                 activebackground='#2980b9',
                                 activeforeground='white')
        sources_browse.grid(row=1, column=2, sticky="w", padx=(10, 0), pady=5)
        
        # Configure grid weights for files frame
        files_frame.grid_columnconfigure(1, weight=1)
        
        # Options section
        options_frame = tk.Frame(main_frame, bg=DnabotApp.colors['background'])
        options_frame.grid(row=5, column=0, columnspan=2, sticky="ew", pady=(0, 20))
        
        # Keep layout option
        self.keep_layout_var = tk.BooleanVar()
        self.keep_layout_var.set(True)
        keep_layout_checkbox = tk.Checkbutton(options_frame, 
                                            text='Keep original CSV layout (preserve empty rows)',
                                            variable=self.keep_layout_var, 
                                            font=DnabotApp.app_font,
                                            fg=DnabotApp.colors['text'],
                                            bg=DnabotApp.colors['background'],
                                            selectcolor=DnabotApp.colors['background'],
                                            activebackground=DnabotApp.colors['background'],
                                            activeforeground=DnabotApp.colors['text'])
        keep_layout_checkbox.pack(anchor="center")
        
        # Buttons section
        buttons_frame = tk.Frame(main_frame, bg=DnabotApp.colors['background'])
        buttons_frame.grid(row=6, column=0, columnspan=2, sticky="ew", pady=(20, 0))
        
        self.quit_status = False
        
        # Style the buttons
        quit_button = tk.Button(buttons_frame, text='QUIT', 
                               fg='white', bg=DnabotApp.colors['accent'],
                               command=self.quitter, 
                               font=DnabotApp.app_font,
                               relief="flat", bd=0,
                               padx=30, pady=10,
                               activebackground='#c0392b',
                               activeforeground='white')
        quit_button.pack(side="left", padx=(0, 10))
        
        generate_button = tk.Button(buttons_frame, text='GENERATE',
                                   fg='white', bg=DnabotApp.colors['success'],
                                   command=self.generate, 
                                   font=DnabotApp.app_font,
                                   relief="flat", bd=0,
                                   padx=30, pady=10,
                                   activebackground='#229954',
                                   activeforeground='white')
        generate_button.pack(side="right")

    def browse_constructs(self):
        """Open file dialog to browse for constructs CSV file"""
        filename = filedialog.askopenfilename(
            title="Select Constructs CSV file",
            filetypes=(("CSV files", "*.csv"), ("All files", "*.*"))
        )
        if filename:
            self.construct_path.set(filename)
    
    def browse_sources(self):
        """Open file dialog to browse for source CSV files"""
        filenames = filedialog.askopenfilenames(
            title="Select Source CSV files",
            filetypes=(("CSV files", "*.csv"), ("All files", "*.*"))
        )
        if filenames:
            # Join multiple file paths with commas
            self.sources_paths.set(','.join(filenames))

    def quitter(self):
        self.quit_status = True
        self.master.quit()

    def generate(self):
        self.etoh_well = self.etoh_well.get()
        self.water_well = self.water_well.get()
        self.thermocycler_gen = self.thermocycler_gen.get()
        self.keep_layout = self.keep_layout_var.get()
        self.construct_path = self.construct_path.get()
        self.sources_paths = self.sources_paths.get()
        self.master.quit()


def main():
    root = tk.Tk()
    dnabotinst = DnabotApp(root)
    root.mainloop()
    root.destroy()
    if dnabotinst.quit_status:
        sys.exit("User specified 'QUIT' during app")
    print('Ethanol well is ', dnabotinst.etoh_well)
    print('Water well is ', dnabotinst.water_well)
    print('Keep layout is ', dnabotinst.keep_layout)
    print('Construct path is ', dnabotinst.construct_path)
    print('Sources paths are ', dnabotinst.sources_paths)


if __name__ == '__main__':
    main()
