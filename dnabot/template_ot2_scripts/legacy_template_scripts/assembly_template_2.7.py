from opentrons import simulate
import numpy as np
# metadata
metadata = {
'protocolName': 'My Protocol',
'author': 'Name <email@address.com>',
'description': 'Simple protocol to get started using OT2',
'apiLevel': '2.7'
}
protocol = simulate.get_protocol_api('2.7')
# protocol run function. the part after the colon lets your editor know
# where to look for autocomplete suggestions
# where to look for autocomplete suggestions

def final_assembly(final_assembly_dict, tiprack_num, tiprack_type="opentrons_96_filtertiprack_10ul"):
            # Constants, we update all the labware name in version 2
            #Tiprack
            CANDIDATE_TIPRACK_SLOTS = ['3', '6', '9', '2', '5', '8', '11']
            PIPETTE_MOUNT = 'right'
            #Plate of sample after  purification
            MAG_PLATE_TYPE = 'nest_96_wellplate_100ul_pcr_full_skirt'
            MAG_PLATE_POSITION = '1'
            #Tuberack
            TUBE_RACK_TYPE = 'opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical'
            TUBE_RACK_POSITION = '7'
            #Destination plate
            DESTINATION_PLATE_TYPE = 'opentrons_96_aluminumblock_generic_pcr_strip_200ul'
            #Temperature control plate
            TEMPDECK_SLOT = '4'
            TEMP = 20
            TOTAL_VOL = 15
            PART_VOL = 1.5
            MIX_SETTINGS = (1, 3)
            tiprack_num=tiprack_num+1
            # Errors
            sample_number = len(final_assembly_dict.keys())
            if sample_number > 96:
                raise ValueError('Final assembly nummber cannot exceed 96.')
                
            slots = CANDIDATE_TIPRACK_SLOTS[:tiprack_num]
            tipracks = [protocol.load_labware(tiprack_type, slot)
                for slot in slots]
            pipette = protocol.load_instrument('p10_single', PIPETTE_MOUNT, tip_racks=tipracks)#old code: pipette = instruments.P10_Single(mount=PIPETTE_MOUNT, tip_racks=tipracks)
            
            # Define Labware and set temperature
            magbead_plate = protocol.load_labware(MAG_PLATE_TYPE, MAG_PLATE_POSITION)
           #old code: magbead_plate = labware.load(MAG_PLATE_TYPE, MAG_PLATE_POSITION)
            tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
           #old code: tube_rack = labware.load(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
            tempdeck = protocol.load_module('tempdeck', TEMPDECK_SLOT)
           #old code: tempdeck = modules.load('tempdeck', TEMPDECK_SLOT)
            destination_plate = tempdeck.load_labware(
            DESTINATION_PLATE_TYPE, TEMPDECK_SLOT)
            tempdeck.set_temperature(TEMP)
           #old code: destination_plate = labware.load(DESTINATION_PLATE_TYPE, TEMPDECK_SLOT, share=True)tempdeck.set_temperature(TEMP)tempdeck.wait_for_temp()
            
             # Master mix transfers
            final_assembly_lens = []
            for values in final_assembly_dict.values():
                final_assembly_lens.append(len(values))
            unique_assemblies_lens = list(set(final_assembly_lens))
            master_mix_well_letters = ['A', 'B', 'C', 'D']
            for x in unique_assemblies_lens:
                master_mix_well = master_mix_well_letters[(x - 1) // 6] + str(x - 1)
                destination_inds = [i for i, lens in enumerate(final_assembly_lens) if lens == x]
                destination_wells = np.array([key for key, value in list(final_assembly_dict.items())])
                destination_wells = list(destination_wells[destination_inds])
                for destination_well in destination_wells:# make tube_rack_wells and destination_plate.wells in the same type
                    pipette.pick_up_tip()
                    pipette.transfer(TOTAL_VOL - x * PART_VOL, tube_rack.wells(master_mix_well),
                                     destination_plate.wells(destination_well), new_tip='never')#transfer water and buffer in the pipette

                    pipette.drop_tip()

            # Part transfers
            for key, values in list(final_assembly_dict.items()):
                for value in values:# magbead_plate.wells and destination_plate.wells in the same type
                    pipette.transfer(PART_VOL, magbead_plate.wells(value),
                                     destination_plate.wells(key), mix_after=MIX_SETTINGS,
                                     new_tip='always')#transfer parts in one tube
            
            tempdeck.deactivate() #stop increasing the temperature
            
final_assembly(final_assembly_dict=final_assembly_dict, tiprack_num=tiprack_num)
            
for line in protocol.commands(): 
    print(line)