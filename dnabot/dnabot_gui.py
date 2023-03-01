# -*- coding: utf-8 -*-
"""
Created on Thu May 30 14:35:26 2019

@author: mh2210, tduigou
"""

from __future__ import annotations  # Enable the "hint" feature for objects

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog


def to_numeric_value(str_value: str):
    if float(str_value).is_integer():
        return int(float(str_value))
    else:
        return float(str_value)


class FileSelector:

    def __init__(self, root, irow, title, multiple_files=False):
        # To be used to get the value from outside
        self.output = None
        self.irow = irow
        # Settings
        self.title = title
        self.multiple_files = multiple_files
        # File path as text
        self.entry_text = tk.Entry(root, state="readonly", width=60)
        self.entry_text.insert(0, '')
        # Browse button
        button = tk.Button(root, text='Browse', command=self.browse_file)
        # Position on grid
        button.grid(row=self.irow, columnspan=2)
        self.irow += 1
        self.entry_text.grid(row=self.irow, columnspan=2)
        self.irow += 1

    def browse_file(self):
        if self.multiple_files:
            output = filedialog.askopenfilenames(  # output is a list
                title=self.title,
                filetypes=(('CSV files', '*.csv'), ('All files', '*.*'))
                )
        else:
            # Output is a str
            output = filedialog.askopenfilename(  # output is a str
                title=self.title,
                filetypes=(('CSV files', '*.csv'), ('All files', '*.*'))
                )
        self.output = output
        self.update_text()
    
    def update_text(self):
        self.entry_text.configure(state='normal')
        self.entry_text.delete(0, "end")
        if isinstance(self.output, tuple):
            self.entry_text.insert(0, ', '.join(self.output))
        else:
            self.entry_text.insert(0, self.output)
        self.entry_text.configure(state='readonly')

    def get(self):
        return self.output


class GUI:

    __APP_TITLE = "DNABot App"
    __APP_FONT = ("Helvetica", 12)
    __TROUGH_WELLS = ['A{}'.format(x + 1) for x in range(12)]

    
    def __init__(
        self,
        root: tk.Tk,
        user_settings = dict
        ) -> GUI:
    
        # The set up the GUI backbone
        self.root = root
        self.canvas = tk.Canvas(self.root, width=650, height=840)
        self.frame = tk.Frame(self.canvas)
        self.vsb = tk.Scrollbar(self.root, orient="vertical", command=self.canvas.yview, width=20)
        self.vsb.pack(side="right", fill="y")
        self.canvas.configure(yscrollcommand=self.vsb.set)
        self.canvas.pack(side="left", fill="both", expand=False)
        self.canvas.create_window((0,0), window=self.frame, anchor="nw", tags="self.frame")

        # 
        self.root.title(GUI.__APP_TITLE)
        self.user_settings = user_settings
        self.quit_status = False

        # Intro ===============================================================
        irow = 0
        intro = tk.Message(
            self.frame,
            text=(
                "Welcome to the dnabot App! Please follow these "
                "instructions to create the 4 DNA-BOT scripts:"),
            width=600)
        intro.grid(row=irow, columnspan=2, padx=5, pady=15)

        # Sep =================================================================
        irow += 1
        self.__add_separator(irow)

        # Ethanol & SOC media =================================================
        irow += 1
        message_1 = tk.Message(
            self.frame,
            text=(
                "1. From the dropdown menus select wells/columns "
                "for ethanol (2_purification script) and SOC media "
                "(4_transformation script)."),
            width=600,
            anchor='w')
        message_1.grid(row=irow, columnspan=2, padx=5, pady=10, sticky='w')

        irow += 1
        etoh_well_label = tk.Label(self.frame, text='Trough well for ethanol during purification:', font=GUI.__APP_FONT)
        etoh_well_label.grid(row=irow, column=0, sticky='e')
        self.etoh_well = tk.StringVar(self.frame)
        self.etoh_well.set(GUI.__TROUGH_WELLS[10])
        etoh_w=tk.OptionMenu(self.frame, self.etoh_well, *tuple(GUI.__TROUGH_WELLS[1:11]))
        etoh_w.grid(row=irow, column=1, sticky=tk.W)
        etoh_w.config(font=GUI.__APP_FONT)

        irow += 1
        soc_column_label = tk.Label(self.frame, text='Deep-well plate column for SOC media during transformation:', font=GUI.__APP_FONT)
        soc_column_label.grid(row=irow, column=0, sticky='e')
        self.soc_column=tk.StringVar(self.frame)
        self.soc_column.set("1")
        soc_w=tk.OptionMenu(self.frame, self.soc_column, *tuple(['{}'.format(x + 1) for x in range(12)]))
        soc_w.grid(row=irow, column=1, sticky='w')
        soc_w.config(font=GUI.__APP_FONT)

        # Sep =================================================================
        irow += 1
        self.__add_separator(irow)

        # Labware IDs =========================================================
        irow += 1
        message_2 = tk.Message(
            self.frame,
            text=(
                "2. Specify the labware IDs to be used. Leave values "
                "as they are to use the default ones."),
            width=600,
            anchor='w')
        message_2.grid(row=irow, columnspan=2, padx=5, pady=10, sticky='w')

        # Opentrons P20 Single-Channel Electronic Pipette
        irow += 1
        self.labware_p10_single_entry = self.__make_labware_entry(
            label="Opentrons P20 Single-Channel Electronic Pipette",
            labware_id='p20_single',
            irow=irow)
        # Opentrons P300 8-Channel Electronic Pipette
        irow += 1
        self.labware_p300_multi_entry = self.__make_labware_entry(
            label="Opentrons P300 8-Channel Electronic Pipette",
            labware_id='p300_multi',
            irow=irow)
        # Opentrons magnetic module
        irow += 1
        self.labware_mag_deck_entry = self.__make_labware_entry(
            label="Opentrons magnetic module",
            labware_id='mag_deck',
            irow=irow)
        # Opentrons 4-in-1 tubes rack for 1.5 ml eppendorf tubes
        irow += 1
        self.labware_24_tuberack_1500ul_entry = self.__make_labware_entry(
            label="Opentrons 4-in-1 tubes rack",
            labware_id='24_tuberack_1500ul',
            irow=irow)
        # Opentrons 10μL tips rack
        irow += 1
        self.labware_96_tiprack_20ul_entry = self.__make_labware_entry(
            label="Opentrons 20μL tips rack",
            labware_id='96_tiprack_20ul',
            irow=irow)
        # Opentrons 300μL tips rack
        irow += 1
        self.labware_96_tiprack_300ul_entry = self.__make_labware_entry(
            label="Opentrons 300μL tips rack",
            labware_id='96_tiprack_300ul',
            irow=irow)

        # Clip reaction source plate (steps: clip)
        irow += 1
        self.labware_clip_source_plate_entry = self.__make_labware_entry(
            label="Clip reaction source plate (steps: clip)",
            labware_id='clip_source_plate',
            irow=irow)

        # Clip reaction plate (steps: clip, purif, assembly)
        irow += 1
        self.labware_clip_plate_entry = self.__make_labware_entry(
            label="Clip reaction plate (steps: clip, purification, assembly)",
            labware_id='clip_plate',
            irow=irow)

        # Mix plate (step: purification)
        irow += 1
        self.labware_mix_plate_entry = self.__make_labware_entry(
            label="Mix plate (step purification)",
            labware_id='mix_plate',
            irow=irow)

        # Final assembly plate (steps: assembly, transformation)
        irow += 1
        self.labware_final_assembly_plate_entry = self.__make_labware_entry(
            label="Final assembly plate (steps: assembly, transformation)",
            labware_id='final_assembly_plate',
            irow=irow)

        # Transformation plate (step: transformation)
        irow += 1
        self.labware_transfo_plate_entry = self.__make_labware_entry(
            label="Transformation plate (step: transformation)",
            labware_id='transfo_plate',
            irow=irow)

        # Transformation plate without thermocycler (step: transformation)
        irow += 1
        self.labware_transfo_plate_wo_thermo_entry = self.__make_labware_entry(
            label="Transformation plate without thermocycler (step: transformation)",
            labware_id='transfo_plate_wo_thermo',
            irow=irow)

        # Agar plate (transformation step)
        irow += 1
        self.agar_plate_entry = self.__make_labware_entry(
            label="Agar plate (transformation step)",
            labware_id='agar_plate',
            irow=irow)

        # Reservoir plate 21 mL 12 channels
        irow += 1
        self.labware_12_reservoir_21000ul_entry = self.__make_labware_entry(
            label="Reservoir plate 21 mL 12 channels",
            labware_id='12_reservoir_21000ul',
            irow=irow)

        # 96 deep well plate 2 mL wells
        irow += 1
        self.labware_96_deepwellplate_2ml_entry = self.__make_labware_entry(
            label="96 deep well plate 2 mL wells",
            labware_id='96_deepwellplate_2ml',
            irow=irow)

        # Corning 12 Well Plate 6.9 mL Flat
        irow += 1
        self.labware_12_corning_wellplate_entry = self.__make_labware_entry(
            label="Corning 12 Well Plate 6.9 mL Flat",
            labware_id="12_corning_wellplate",
            irow=irow)
        
        # Sep =================================================================
        irow += 1
        self.__add_separator(irow)

        # Parameters for the clip reaction step ===============================
        irow += 1
        message_3 = tk.Message(
            self.frame,
            text="3. Specify parameters for the clip reaction step.",
            width=600,
            anchor='w')
        message_3.grid(row=irow, columnspan=2, padx=5, pady=10, sticky='w')
        irow += 1
        self.param_clip_thermo_lid_closed = self.__make_parameter_entry(
            label="How long to keep at 4°C samples at the end of the execution? (in minutes)",
            parameter_id="clip_keep_sample_overnight",
            irow=irow)

        # Sep =================================================================
        irow += 1
        self.__add_separator(irow)

        # Parameters for the purification step ================================
        irow += 1
        message_4 = tk.Message(
            self.frame,
            text="3. Specify parameters for the purification step.",
            width=600,
            anchor='w')
        message_4.grid(row=irow, columnspan=2, padx=5, pady=10, sticky='w')
        irow += 1
        self.param_purif_magdeck_height = self.__make_parameter_entry(
            label="Magnetic module height (mm)",
            parameter_id="purif_magdeck_height",
            irow=irow)
        irow += 1
        self.param_purif_wash_time = self.__make_parameter_entry(
            label="Washing time (min)",
            parameter_id="purif_wash_time",
            irow=irow)
        irow += 1
        self.param_purif_bead_ratio = self.__make_parameter_entry(
            label="Bead ratio",
            parameter_id="purif_bead_ratio",
            irow=irow)
        irow += 1
        self.param_purif_incubation_time = self.__make_parameter_entry(
            label="Incubation time (min)",
            parameter_id="purif_incubation_time",
            irow=irow)
        irow += 1
        self.param_purif_settling_time = self.__make_parameter_entry(
            label="Settling time (min)",
            parameter_id="purif_settling_time",
            irow=irow)
        irow += 1
        self.param_purif_drying_time = self.__make_parameter_entry(
            label="Drying time (min)",
            parameter_id="purif_drying_time",
            irow=irow)
        irow += 1
        self.param_purif_elution_time = self.__make_parameter_entry(
            label="Elution time (min)",
            parameter_id="purif_elution_time",
            irow=irow)

        # Sep =================================================================
        irow += 1
        self.__add_separator(irow)

        # Parameters for the transformation step ==============================
        irow += 1
        message_5 = tk.Message(
            self.frame,
            text="4. Specify parameters for the transformation step.",
            width=600,
            anchor='w')
        message_5.grid(row=irow, columnspan=2, padx=5, pady=10, sticky='w')
        irow += 1
        self.param_transfo_incubation_temp = self.__make_parameter_entry(
            label="Incubation temperature (°C)",
            parameter_id="transfo_incubation_temp",
            irow=irow)
        irow += 1
        self.param_transfo_incubation_time = self.__make_parameter_entry(
            label="Incubation time (min)",
            parameter_id="transfo_incubation_time",
            irow=irow)

        # Sep =================================================================
        irow += 1
        self.__add_separator(irow)

        # Construct CSV file ==================================================
        irow += 1
        message_6 = tk.Message(
            self.frame,
            text="5. Select the CSV file describing constructs.",
            width=600,
            anchor='w')
        message_6.grid(row=irow, columnspan=2, padx=5, pady=10, sticky='w')
        irow += 1
        self.construct_file_selector = FileSelector(self.frame, irow, title='Construct CSV file', multiple_files=False)
        irow = self.construct_file_selector.irow

        # Sep =================================================================
        irow += 1
        self.__add_separator(irow)

        # Source CSV files ====================================================
        irow += 1
        message_7 = tk.Message(
            self.frame,
            text=(
                "6. Select up to 6 csv files describing plates "
                "containing BASIC parts and linkers. If all files "
                "are not within one folder, absolute paths should "
                "be given."),
            width=600,
            anchor='w')
        message_7.grid(row=irow, columnspan=2, padx=5, pady=10, sticky='w')
        irow += 1
        self.source_files_selector = FileSelector(self.frame, irow, title='Source CSV files', multiple_files=True)
        irow = self.source_files_selector.irow

        # Sep =================================================================
        irow += 1
        self.__add_separator(irow)

        # # White space
        # irow += 1
        # spacer = tk.Label(self.frame, text="", font=GUI.__APP_FONT)
        # spacer.grid(row=irow, columnspan=2, padx=5, pady=10)

        # Quit and generate buttons
        irow += 1
        quit_button = tk.Button(self.frame, text='QUIT', fg='red', command=self.quit, font=GUI.__APP_FONT)
        quit_button.grid(row=irow, column=0, pady=10)
        generate_button=tk.Button(self.frame, text='GENERATE', command=self.generate, font=GUI.__APP_FONT)
        generate_button.grid(row=irow, column=1, pady=10)

        # Upate scroll region info
        self.root.update()  # Required to refresh canvas scrollreion
        self.canvas.configure(scrollregion = self.canvas.bbox("all"))

        # 
        self.frame.mainloop()


    def quit(self):
        self.quit_status=True
        self.root.quit()

    def generate(self):
        # Step 1
        self.user_settings['etoh_well'] = self.etoh_well.get()
        self.user_settings['soc_column'] = self.soc_column.get()
        # Step 2
        self.user_settings['labwares']['p20_single']['id'] = self.labware_p10_single_entry.get()
        self.user_settings['labwares']['p300_multi']['id'] = self.labware_p300_multi_entry.get()
        self.user_settings['labwares']['mag_deck']['id'] = self.labware_mag_deck_entry.get()
        self.user_settings['labwares']['24_tuberack_1500ul']['id'] = self.labware_24_tuberack_1500ul_entry.get()
        self.user_settings['labwares']['96_tiprack_20ul']['id'] = self.labware_96_tiprack_20ul_entry.get()
        self.user_settings['labwares']['96_tiprack_300ul']['id'] = self.labware_96_tiprack_300ul_entry.get()

        self.user_settings['labwares']['clip_source_plate']['id'] = self.labware_clip_source_plate_entry.get()
        self.user_settings['labwares']['clip_plate']['id'] = self.labware_clip_plate_entry.get()
        self.user_settings['labwares']['mix_plate']['id'] = self.labware_mix_plate_entry.get()
        self.user_settings['labwares']['final_assembly_plate']['id'] = self.labware_final_assembly_plate_entry.get()
        self.user_settings['labwares']['transfo_plate']['id'] = self.labware_transfo_plate_entry.get()
        self.user_settings['labwares']['transfo_plate_wo_thermo']['id'] = self.labware_transfo_plate_wo_thermo_entry.get()

        self.user_settings['labwares']['agar_plate']['id'] = self.agar_plate_entry.get()
        self.user_settings['labwares']['12_reservoir_21000ul']['id'] = self.labware_12_reservoir_21000ul_entry.get()
        self.user_settings['labwares']['96_deepwellplate_2ml']['id'] = self.labware_96_deepwellplate_2ml_entry.get()
        self.user_settings['labwares']['12_corning_wellplate']['id'] = self.labware_12_corning_wellplate_entry.get()
        # Step 3
        self.user_settings['parameters']['purif_magdeck_height']['value'] = to_numeric_value(self.param_purif_magdeck_height.get())
        self.user_settings['parameters']['purif_wash_time']['value'] = to_numeric_value(self.param_purif_wash_time.get())
        self.user_settings['parameters']['purif_bead_ratio']['value'] = to_numeric_value(self.param_purif_bead_ratio.get())
        self.user_settings['parameters']['purif_incubation_time']['value'] = to_numeric_value(self.param_purif_incubation_time.get())
        self.user_settings['parameters']['purif_settling_time']['value'] = to_numeric_value(self.param_purif_settling_time.get())
        self.user_settings['parameters']['purif_drying_time']['value'] = to_numeric_value(self.param_purif_drying_time.get())
        self.user_settings['parameters']['purif_elution_time']['value'] = to_numeric_value(self.param_purif_elution_time.get())
        # Step 4
        self.user_settings['parameters']['transfo_incubation_temp']['value'] = to_numeric_value(self.param_transfo_incubation_temp.get())
        self.user_settings['parameters']['transfo_incubation_time']['value'] = to_numeric_value(self.param_transfo_incubation_time.get())
        # Step 5
        self.user_settings['construct_path'] = self.construct_file_selector.get()
        # Step 6
        self.user_settings['sources_paths'] = self.source_files_selector.get()
        self.root.quit()

    def __make_labware_entry(self, label, labware_id, irow):
        labware_label = tk.Label(self.frame, text=label, font=GUI.__APP_FONT)
        labware_label.grid(row=irow, column=0, sticky='e')
        labware_entry = tk.Entry(self.frame, width=28)
        labware_entry.insert(0, self.user_settings["labwares"][labware_id]['id'])
        labware_entry.grid(row=irow, column=1, sticky='w')
        return labware_entry

    def __make_parameter_entry(self, label, parameter_id, irow, parameter_value="value"):
        parameter_label = tk.Label(self.frame, text=label, font=GUI.__APP_FONT)
        parameter_label.grid(row=irow, column=0, sticky='e')
        parameter_entry = tk.Entry(self.frame, width=28)
        parameter_entry.insert(0, self.user_settings["parameters"][parameter_id][parameter_value])
        parameter_entry.grid(row=irow, column=1, sticky='w')
        return parameter_entry

    def __add_separator(self, irow):
        ttk.Separator(
            self.frame,
            orient=tk.HORIZONTAL
        ).grid(
            row=irow,
            columnspan=2,
            sticky="ew"
        )