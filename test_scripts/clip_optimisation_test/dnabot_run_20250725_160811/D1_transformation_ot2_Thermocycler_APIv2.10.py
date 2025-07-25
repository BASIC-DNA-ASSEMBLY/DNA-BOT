from opentrons import protocol_api
import json
import time

metadata = {
    'apiLevel': '2.10',
    'protocolName': 'Transformation Protocol v2',
    'description': 'Enhanced transformation protocol using thermocycler module',
    'author': 'Liam Hallett'
}

# Load transformation data from JSON file
# This will be replaced by the parser with embedded JSON data
transformation_dict = {
    "transformation_number": 16,
    "source_wells": ["A2", "B2", "C2", "D2", "E2", "F2", "G2", "H2", "A3", "B3", "C3", "D3", "E3", "F3", "G3", "H3"],
    "source_plates": ["1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1"],
    "destination_wells": ["A2", "B2", "C2", "D2", "E2", "F2", "G2", "H2", "A3", "B3", "C3", "D3", "E3", "F3", "G3", "H3"],
    "dna_volumes": [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    "cell_volumes": [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30],
    "soc_volumes": [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
    "plating_volumes": [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]
}

# Thermocycler generation setting
# This will be replaced by the parser with embedded thermocycler generation
thermocycler_gen = 'gen2'

# opentrons_simulate.exe dnabot\template_ot2_scripts\transformation_template_TC_APIv2.10.py --custom-labware-path 'labware\Labware definitions'

# Example dictionary structure for transformation data:
# transformation_dict = {
#     'transformation_number': 24,
#     'source_wells': ['A1', 'A2', 'A3', 'B1', 'B2', 'B3', ...],
#     'source_plates': ['1', '1', '1', '1', '1', '1', ...],
#     'destination_wells': ['A1', 'A2', 'A3', 'B1', 'B2', 'B3', ...],
#     'dna_volumes': [3, 3, 3, 3, 3, 3, ...],
#     'cell_volumes': [30, 30, 30, 30, 30, 30, ...],
#     'soc_volumes': [100, 100, 100, 100, 100, 100, ...],
#     'plating_volumes': [100, 100, 100, 100, 100, 100, ...]
# }

def run(protocol: protocol_api.ProtocolContext):
    """Main transformation protocol function."""
    
    # Extract data from embedded dictionary
    transformation_number = transformation_dict['transformation_number']
    source_wells = transformation_dict['source_wells']
    source_plates = transformation_dict['source_plates']
    destination_wells = transformation_dict['destination_wells']
    dna_volumes = transformation_dict['dna_volumes']
    cell_volumes = transformation_dict['cell_volumes']
    soc_volumes = transformation_dict['soc_volumes']
    plating_volumes = transformation_dict['plating_volumes']
    
    # Load hardware
    left_pipette = protocol.load_instrument('p300_multi_gen2', 'left')
    right_pipette = protocol.load_instrument('p20_single_gen2', 'right')
    
    # Load tipracks
    tiprack300_1 = protocol.load_labware('opentrons_96_tiprack_300ul', 6)
    tiprack300_1_count = 0
    tiprack20_1 = protocol.load_labware('opentrons_96_tiprack_20ul', 9)
    tiprack20_1_count = 0
    
    # Load plates and vessels
    source_plate_a = protocol.load_labware('biorad_96_wellplate_200ul_pcr', 1)
    source_plate_b = protocol.load_labware('biorad_96_wellplate_200ul_pcr', 2)
    plate_2 = protocol.load_labware('corning_12_wellplate_6.9ml_flat', 5)
    tube_rack = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', 4)
    
    # Load thermocycler module
    if thermocycler_gen == 'gen1':
        thermocycler = protocol.load_module('thermocycler', 7)
    else:  # gen2
        thermocycler = protocol.load_module('thermocyclerModuleV2', 7)
    TC_plate = thermocycler.load_labware('biorad_96_wellplate_200ul_pcr')
    
    # Helper functions
    def flash(n):
        """Flashes light n times to prompt user input."""
        for i in range(n):
            protocol.set_rail_lights(False)
            time.sleep(0.13)
            protocol.set_rail_lights(True)
            time.sleep(0.13)
    
    def iterate_tip_count(tip_count, pipette, tiprack):
        """Iterates tip count and prompts for new tip box if needed."""
        pipette.drop_tip()
        tip_box_capacity = 95  # 96-1 due to zero indexing
        
        if tip_count == tip_box_capacity:
            comment = "Change tip box"
            protocol.pause(comment)
            flash(5)
            return 0
        else:
            return tip_count + 1
    
    def multi2single_tip_counter(n):
        """Returns tip position for multi-channel single tip usage."""
        row_dict = {0: "A", 1: "B", 2: "C", 3: "D", 4: "E", 5: "F", 6: "G", 7: "H"}
        row = row_dict[7 - n % 8]
        col = 1 + (n // 8)
        return row + str(col)
    
    def tube_counter(n):
        """Returns tube position in tube rack."""
        row_dict = {0: "A", 1: "B", 2: "C", 3: "D", 4: "E", 5: "F"}
        row = row_dict[n % 4]
        col = 1 + n // 4
        return row + str(col)
    
    def plate_counter(n):
        """Returns plate position for plating."""
        row_dict = {0: "A", 1: "B", 2: "C"}
        row = row_dict[n % 3]
        col = 1 + n // 3
        return row + str(col)
    
    # Calculate requirements
    cell_tube_vol = 1000  # volume of cells per aliquot
    cell_tube_capacity = cell_tube_vol // 30  # transfers per tube (assuming 30µl per transformation)
    cell_req_tubes = (transformation_number - 1) // cell_tube_capacity + 1
    
    SOC_tube_vol = 1500  # vol per tube of SOC
    SOC_tube_capacity = SOC_tube_vol // 100  # transfers per tube (assuming 100µl per transformation)
    req_SOC_tubes = (transformation_number - 1) // SOC_tube_capacity + 1
    
    plate_capacity = 12  # number of wells available on plate
    req_plates = (transformation_number - 1) // plate_capacity + 1
    
    # Protocol execution
    protocol.set_rail_lights(True)
    thermocycler.set_block_temperature(4)
    
    comment = f"Required:\n {cell_req_tubes} tube(s) of competent cells,\n{req_SOC_tubes} tube(s) of SOC,\n{req_plates} plate(s)"
    protocol.pause(comment)
    flash(5)
    
    # Dispense competent cells
    protocol.comment("Dispensing competent cells")
    for i in range(transformation_number):
        tip_count = i % 10  # 300µl pipette can do ~10 transfers of 30µl
        
        if tip_count == 0:
            left_pipette.pick_up_tip(tiprack300_1[multi2single_tip_counter(tiprack300_1_count)])
        
        cell_tube_count = i // cell_tube_capacity
        left_pipette.aspirate(cell_volumes[i], tube_rack[tube_counter(cell_tube_count)])
        
        dispense_count = max(i, 1) % 10
        if dispense_count == 9 or i == transformation_number - 1:
            transfers_completed = (i // 10) * 10
            req_transfers = i - transfers_completed + 1
            
            for transfer in range(req_transfers):
                left_pipette.dispense(cell_volumes[transfers_completed + transfer], 
                                    TC_plate[destination_wells[transfers_completed + transfer]])
            
            left_pipette.touch_tip()
            tiprack300_1_count = iterate_tip_count(tiprack300_1_count, left_pipette, tiprack300_1)
    
    # Transfer DNA
    protocol.comment("Transferring DNA")
    for i in range(transformation_number):
        right_pipette.pick_up_tip(tiprack20_1[f"A{tiprack20_1_count + 1}"])
        
        # Determine source plate
        source_plate = source_plate_a if source_plates[i] == '1' else source_plate_b
        
        right_pipette.aspirate(dna_volumes[i], source_plate[source_wells[i]])
        right_pipette.dispense(dna_volumes[i], TC_plate[destination_wells[i]])
        right_pipette.mix(dna_volumes[i], 10, TC_plate[destination_wells[i]])
        right_pipette.touch_tip()
        
        tiprack20_1_count = iterate_tip_count(tiprack20_1_count, right_pipette, tiprack20_1)
    
    # Heat shock protocol
    protocol.comment("Starting heat shock protocol")
    thermocycler.close_lid()
    thermocycler.set_block_temperature(4, hold_time_minutes=30, block_max_volume=55)
    thermocycler.set_block_temperature(42, hold_time_seconds=30, block_max_volume=55)
    thermocycler.set_block_temperature(4, hold_time_seconds=90, block_max_volume=55)
    
    # Add SOC
    protocol.comment("Adding SOC")
    thermocycler.open_lid()
    thermocycler.set_block_temperature(37)
    
    SOC_tube_count = 4  # starts at A2 (tube 4)
    for i in range(transformation_number):
        step = max(i, 1) % SOC_tube_capacity
        
        left_pipette.pick_up_tip(tiprack300_1[multi2single_tip_counter(tiprack300_1_count)])
        
        if step == 0:
            SOC_tube_count += 1
        
        left_pipette.aspirate(soc_volumes[i], tube_rack[tube_counter(SOC_tube_count)])
        left_pipette.dispense(soc_volumes[i], TC_plate[destination_wells[i]])
        left_pipette.mix(2, 50, TC_plate[destination_wells[i]])
        
        tiprack300_1_count = iterate_tip_count(tiprack300_1_count, left_pipette, tiprack300_1)
    
    # Recovery incubation
    protocol.comment("Recovery incubation")
    thermocycler.close_lid()
    thermocycler.set_block_temperature(37, hold_time_minutes=60, block_max_volume=150)
    
    # Plating
    protocol.comment("Plating transformations")
    thermocycler.open_lid()
    thermocycler.set_block_temperature(37)
    
    for i in range(transformation_number):
        step = i % plate_capacity
        
        if step == 0:
            plate = str(i // plate_capacity + 1)
            comment = f"Place plate {plate} onto robot at deck position 5"
            protocol.pause(comment)
            flash(5)
        
        left_pipette.pick_up_tip(tiprack300_1[multi2single_tip_counter(tiprack300_1_count)])
        left_pipette.aspirate(plating_volumes[i], TC_plate[destination_wells[i]])
        left_pipette.blow_out(plate_2[plate_counter(step)].top())
        left_pipette.touch_tip()
        
        tiprack300_1_count = iterate_tip_count(tiprack300_1_count, left_pipette, tiprack300_1)
    
    # Finish protocol
    thermocycler.deactivate()
    protocol.set_rail_lights(False)
    protocol.comment("Transformation protocol complete")
