# -*- coding: utf-8 -*-
"""
Created on Thu May 30 14:35:26 2019

@author: mh2210
"""

from __future__ import annotations  # Enable the "hint" feature for objects

import tkinter as tk
from tkinter import filedialog


def to_numeric_value(str_value: str):
    if int(str_value) == float(str_value):
        return int(str_value)
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

    __APP_TITLE = "dnabot app"
    __INTRO_TEXT = ("Welcome to the dnabot App! Please follow these "
                    "instructions to create the 4 DNA-BOT scripts:"
                   )
    __INSTRUCTION_STEP_1 = ("1. From the dropdown menus select wells/columns "
                            "for ethanol (2_purification script) and SOC media "
                            "(4_transformation script)."
                           )
    __INSTRUCTION_STEP_2 = ("2. Specify the labware IDs to be used. Leave values "
                            "as they are to use the default ones.")
    __INSTRUCTION_STEP_3 = ("3. Specify parameters for the purification step.")
    __INSTRUCTION_STEP_4 = ("4. Specify parameters for the transformation step.")
    __INSTRUCTION_STEP_5 = ("5. Select the CSV file describing constructs.")
    __INSTRUCTION_STEP_6 = ("6. Select up to 6 csv files describing plates "
                            "containing BASIC parts and linkers. If all files "
                            "are not within one folder, absolute paths should "
                            "be given.")
    __APP_FONT = ("Helvetica", 12)
    __TROUGH_WELLS = ['A{}'.format(x + 1) for x in range(12)]

    
    def __init__(
        self,
        root: tk.Tk,
        user_settings = dict
        ) -> GUI:
     
        self.root = root
        self.root.lift()
        self.root.title(GUI.__APP_TITLE)
        self.user_settings = user_settings
        self.quit_status = False

        # Intro
        irow = 0
        intro = tk.Message(self.root, text=GUI.__INTRO_TEXT, width=600)
        intro.grid(row=irow, columnspan=2, padx=5, pady=15)

        # Step 1 -- Ethanol & SOC media
        irow += 1
        step_1 = tk.Message(self.root, text=GUI.__INSTRUCTION_STEP_1, width=600, anchor='w')
        step_1.grid(row=irow, columnspan=2, padx=5, pady=10, sticky='w')

        irow += 1
        etoh_well_label = tk.Label(self.root, text='Trough well for ethanol during purification:', font=GUI.__APP_FONT)
        etoh_well_label.grid(row=irow, column=0, sticky='e')
        self.etoh_well = tk.StringVar(root)
        self.etoh_well.set(GUI.__TROUGH_WELLS[10])
        etoh_w=tk.OptionMenu(root, self.etoh_well, *tuple(GUI.__TROUGH_WELLS[1:11]))
        etoh_w.grid(row=irow, column=1, sticky=tk.W)
        etoh_w.config(font=GUI.__APP_FONT)

        irow += 1
        soc_column_label = tk.Label(self.root, text='Deep-well plate column for SOC media during transformation:', font=GUI.__APP_FONT)
        soc_column_label.grid(row=irow, column=0, sticky='e')
        self.soc_column=tk.StringVar(root)
        self.soc_column.set("1")
        soc_w=tk.OptionMenu(root, self.soc_column, *tuple(['{}'.format(x + 1) for x in range(12)]))
        soc_w.grid(row=irow, column=1, sticky='w')
        soc_w.config(font=GUI.__APP_FONT)

        # Step 2 -- Labware IDs
        irow += 1
        step_2 = tk.Message(self.root, text=GUI.__INSTRUCTION_STEP_2, width=600, anchor='w')
        step_2.grid(row=irow, columnspan=2, padx=5, pady=10, sticky='w')

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

        # 96 well rigid PCR plate (clip and transformation steps)
        irow += 1
        self.labware_96_wellplate_200ul_pcr_step_14_entry = self.__make_labware_entry(
            label="96 well rigid PCR plate (clip and transformation steps)",
            labware_id='96_wellplate_200ul_pcr_step_14',
            irow=irow)

        # 96 well rigid PCR plate (purification step)
        irow += 1
        self.labware_96_wellplate_200ul_pcr_step_23_entry = self.__make_labware_entry(
            label="96 well rigid PCR plate (purification and assembly steps)",
            labware_id='96_wellplate_200ul_pcr_step_23',
            irow=irow)

        # Agar plate (transformation step)
        irow += 1
        self.agar_plate_step_4_entry = self.__make_labware_entry(
            label="Agar plate (transformation step)",
            labware_id='agar_plate_step_4',
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
        
        # Step 3 -- Parameters for the purification step
        irow += 1
        step_3 = tk.Message(self.root, text=GUI.__INSTRUCTION_STEP_3, width=600, anchor='w')
        step_3.grid(row=irow, columnspan=2, padx=5, pady=10, sticky='w')
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

        # Step 4 -- Parameters for the transformation step
        irow += 1
        step_4 = tk.Message(self.root, text=GUI.__INSTRUCTION_STEP_4, width=600, anchor='w')
        step_4.grid(row=irow, columnspan=2, padx=5, pady=10, sticky='w')
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

        # Step 5 -- Construct CSV file
        irow += 1
        step_5 = tk.Message(self.root, text=GUI.__INSTRUCTION_STEP_5, width=600, anchor='w')
        step_5.grid(row=irow, columnspan=2, padx=5, pady=10, sticky='w')
        irow += 1
        self.construct_file_selector = FileSelector(root, irow, title='Construct CSV file', multiple_files=False)
        irow = self.construct_file_selector.irow

        # Step 6 -- Source CSV files
        irow += 1
        step_6 = tk.Message(self.root, text=GUI.__INSTRUCTION_STEP_6, width=600, anchor='w')
        step_6.grid(row=irow, columnspan=2, padx=5, pady=10, sticky='w')
        irow += 1
        self.source_files_selector = FileSelector(root, irow, title='Source CSV files', multiple_files=True)
        irow = self.source_files_selector.irow

        # White space
        irow += 1
        spacer = tk.Label(root, text="", font=GUI.__APP_FONT)
        spacer.grid(row=irow, columnspan=2, padx=5, pady=10)

        # Quit and generate buttons
        irow += 1
        quit_button = tk.Button(self.root, text='QUIT', fg='red', command=self.quit, font=GUI.__APP_FONT)
        quit_button.grid(row=irow, column=0, pady=10)
        generate_button=tk.Button(root, text='GENERATE', command=self.generate, font=GUI.__APP_FONT)
        generate_button.grid(row=irow, column=1, pady=10)

        root.mainloop()


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
        self.user_settings['labwares']['96_wellplate_200ul_pcr_step_14']['id'] = self.labware_96_wellplate_200ul_pcr_step_14_entry.get()
        self.user_settings['labwares']['96_wellplate_200ul_pcr_step_23']['id'] = self.labware_96_wellplate_200ul_pcr_step_23_entry.get()
        self.user_settings['labwares']['agar_plate_step_4']['id'] = self.agar_plate_step_4_entry.get()
        self.user_settings['labwares']['12_reservoir_21000ul']['id'] = self.labware_12_reservoir_21000ul_entry.get()
        self.user_settings['labwares']['96_deepwellplate_2ml']['id'] = self.labware_96_deepwellplate_2ml_entry.get()
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
        labware_label = tk.Label(self.root, text=label, font=GUI.__APP_FONT)
        labware_label.grid(row=irow, column=0, sticky='e')
        labware_entry = tk.Entry(self.root, width=30)
        labware_entry.insert(0, self.user_settings['labwares'][labware_id]['id'])
        labware_entry.grid(row=irow, column=1, sticky='w')
        return labware_entry

    def __make_parameter_entry(self, label, parameter_id, irow):
        parameter_label = tk.Label(self.root, text=label, font=GUI.__APP_FONT)
        parameter_label.grid(row=irow, column=0, sticky='e')
        parameter_entry = tk.Entry(self.root, width=30)
        parameter_entry.insert(0, self.user_settings['parameters'][parameter_id]['value'])
        parameter_entry.grid(row=irow, column=1, sticky='w')
        return parameter_entry


# def main():
#     #    root = tk.Tk()
#     #    construct_path = UserDefinedPaths(root, 'Construct csv file')
#     #    root.destroy()
#     #    constructs_list = basic_ot2_generator_v11_beta.generate_constructs_list(
#     #            construct_path.output)
#     #    print(constructs_list)

#     default_user_settings = {"labwares": {"p10_single": {"id": "p20_single_gen2"}, "p300_multi": {"id": "p300_multi_gen2"}, "mag_deck": {"id": "magdeck"}, "24_tuberack_1500ul": {"id": "e14151500starlab_24_tuberack_1500ul"}, "96_tiprack_20ul": {"id": "opentrons_96_tiprack_20ul"}, "96_tiprack_300ul": {"id": "opentrons_96_tiprack_300ul"}, "96_wellplate_200ul_pcr": {"id": "4ti0960rig_96_wellplate_200ul"}, "12_reservoir_21000ul": {"id": "4ti0131_12_reservoir_21000ul"}, "96_deepwellplate_2ml": {"id": "4ti0136_96_wellplate_2200ul"}}}

#     root = tk.Tk()
#     dnabotinst = GUI(root, default_user_settings)
#     root.mainloop()
#     root.destroy()
#     if dnabotinst.quit_status:
#         sys.exit("User specified 'QUIT' during app")
#     print('Ethanol well is ', dnabotinst.etoh_well)
#     print('SOC column is ', dnabotinst.soc_column)


# if __name__ == '__main__':
#     main()
