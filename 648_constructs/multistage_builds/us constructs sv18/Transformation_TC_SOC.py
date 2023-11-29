from opentrons import protocol_api, simulate
import opentrons, time

metadata = {'apiLevel': '2.10'}

def counter(rows):
    def inner(n):
        """
        Takes a number 'n' between 0 and 95, representing the nth tip in the box.
        Returns a well in the form 'A1'.
        """
        row_dict = {0: "A", 1: "B", 2: "C", 3: "D", 4: "E", 5: "F", 6: "G", 7: "H"}
        row = row_dict[n % rows]
        col = 1 + n // rows
        return row + str(col)
    return inner

tip_counter = counter(8)
tube_counter = counter(4)
plate_counter = counter(3)

def counter_T(cols):
    def inner(n):
        """
        Takes a number 'n' between 0 and 95, representing the nth tip in the box.
        Returns a well in the form 'A1'. 
        Tip order is transposed (i.e. A1, A2, ...)
        """
        row_dict = {0: "A", 1: "B", 2: "C", 3: "D", 4: "E", 5: "F", 6: "G", 7: "H"}
        row = row_dict[n // cols]
        col = 1 + n % cols
        return row + str(col)
    return inner

tip_counter_T = counter_T(18)
plate_counter_T = counter_T(4)

def multi2single_tip_counter(n):
    """
    Takes a number 'n' between 0 and 95, representing the nth tip in the box.
    Returns a well in the form 'A1'. The tips are taken in an order that 
    allows the mutli channel to pick up one tip at a time. This is so that only
    one channel of the multi is used limiting the effects of variations between
    channels.

    NB: DO NOT put the tip box or destination plates on positions 1-3 due to risk of collision.
    """
    row_dict = {0: "A", 1: "B", 2: "C", 3: "D", 4: "E", 5: "F", 6: "G", 7: "H"}
    row = row_dict[7 - n % 8]
    col = 1 + (n // 8)
    return row + str(col)
    



# t_array_a = ["A1", "A2", "A3"]#, "B7", "B8", "C1", "E10", "E12"]
t_array_a = [tip_counter(i) for i in range(5)]         # source a wells list
t_array_b = [tip_counter(i) for i in range(0)]          # source b wells list for second plate if needed (deck space 2)

# print([tip_counter(i) for i in range(25)])
print(t_array_a)



def run(protocol: protocol_api.ProtocolContext, t_array_a = t_array_a, t_array_b = t_array_b):

    ###################
    #### Functions ####
    ###################

    def flash(n):
        ''' Flashes light three times. Intended to prompt user input. '''
        for i in range(n):
            protocol.set_rail_lights(False)
            time.sleep(0.13)
            protocol.set_rail_lights(True)
            time.sleep(0.13)

    def iterate_tip_count(tip_count, pipette):
        ''' Iterates tip count by 1 whilst also checking whether tip box capacity (i.e. 96) has been reached. 
            If so, user is prompted for new box. This approach is a workaround for tip box limitations'''

        pipette.drop_tip()

        tip_box_capacity = 95 # 96-1 due to zero indexing

        if tip_count == tip_box_capacity:
            comment = "Change tip box"
            protocol.pause(comment)
            flash(5)

            return 0

        else:
            return tip_count + 1

    ##################
    #### Hardware ####
    ##################

    # Pipettes
    left = protocol.load_instrument('p300_multi_gen2', 'left')
    right = protocol.load_instrument('p20_single_gen2', 'right')

    # Tipracks
    tiprack300_1 = protocol.load_labware('opentrons_96_tiprack_300ul', 6)
    tiprack300_1_count = 0

    tiprack20_1 = protocol.load_labware('opentrons_96_tiprack_20ul', 9)
    tiprack20_1_count = 0
    
    # Plates and vessels 
    source_plate_a = protocol.load_labware('biorad_96_wellplate_200ul_pcr', 1)
    source_plate_b = protocol.load_labware('biorad_96_wellplate_200ul_pcr', 2)
    plate_2 = protocol.load_labware('corning_12_wellplate_6.9ml_flat', 5)
    tube_rack = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', 4)

    # Thermocycler Module GEN1
    thermocycler = protocol.load_module('thermocycler', 7)
    TC_plate = thermocycler.load_labware('biorad_96_wellplate_200ul_pcr')

    # Temperature module
    # temperature_module = protocol.load_module('temperature module gen2', 4)
    # tube_rack = temperature_module.load_labware('opentrons_24_aluminumblock_generic_2ml_screwcap')

    ###################
    #### Constants ####
    ###################

    # Dealing with multiple plates
    total_len = len(t_array_a) + len(t_array_b)
    t_array = [tip_counter(i) for i in range(total_len)]    # destination wells list

    source_dict = {}
    well_dict   = {}

    for i in range(len(t_array)):
        source_dict[i]  = (source_plate_b, source_plate_a)[i < len(t_array_a)]
        well_dict[i]    = (tip_counter(i-len(t_array_a)), tip_counter(i))[i < len(t_array_a)]


    # DNA
    DNA_vol = 3                                                     # vol of DNA per transformation

    # Cells
    cells_vol = 30                                                  # vol of cells per transformation
    p300_vol = 300                                                  # pipette volume
    cell_tube_vol = 1000                                             # volume of cells per aliquot
    air_gap = 0

    cell_tip_capacity = p300_vol // (cells_vol + air_gap)           # tip multi_transfer capacity (ie number of transfers possible with one tip)
    # cell_req_tips = (len(t_array) // cell_tip_capacity) + 1         # tips required to complete transfer of cells

    cell_tube_capacity = cell_tube_vol // cells_vol                 # transfers per tube containing 750ul of cells
    cell_req_tubes = (len(t_array) - 1)  // cell_tube_capacity + 1  # tubes required for total volume of cells

    # SOC
    SOC_transfer_vol = 100                                          # vol of SOC to add per transformation
    SOC_tube_vol = 1500                                             # vol per tube of SOC
    SOC_tube_count = 4                                              # current SOC tube counter (starts at 4 - A2)

    SOC_tube_capacity = SOC_tube_vol // SOC_transfer_vol            # transfers per tube of SOC    
    req_SOC_tubes = (len(t_array) - 1) // SOC_tube_capacity + 1     # tubes required for total volume of SOC
    # SOC_req_tips = len(t_array)                                     # tips required to complete transfer of SOC

    # Plating
    plating_vol = 100                                               # volume of transformed mix to plate out
    plate_capacity = 12                                             # number of wells available on plate
    req_plates = (len(t_array) - 1) // plate_capacity + 1           # plates required

    comment = "Required:\n " + str(cell_req_tubes) + " tube(s) of competent cells,\n" + str(req_SOC_tubes) + " tube(s) of SOC,\n" + str(req_plates) + " plate(s),\n"
    protocol.pause(comment)
    flash(5)

    ##################
    #### Protocol ####
    ##################

    protocol.set_rail_lights(True)
    thermocycler.set_block_temperature(4)
    # temperature_module.set_temperature(4)

    # # Dispense 30ul of competent cells
    # for transformant in range(len(t_array)):                                # NB: uses one tip for multiple transfers
        
    #     tip_count = transformant % cell_tip_capacity                        # step of the multiple transfer tip capacity

    #     if tip_count == 0:                                                  # triggers new tip, includes initial tip
    #         left.pick_up_tip(tiprack300_1[multi2single_tip_counter(tiprack300_1_count)])
        
    #     cell_tube_count = transformant // cell_tube_capacity

    #     # left.air_gap(air_gap/2)
    #     left.aspirate(cells_vol, tube_rack[tube_counter(cell_tube_count)])  # aspirate cells with air gap
    #     # left.air_gap(air_gap/2)

    #     dispense_count = max(transformant, 1) % cell_tip_capacity           # finds the nth transformant, not including the first, triggers dispensing when == 0 

    #     if dispense_count == cell_tip_capacity - 1 or transformant == len(t_array) - 1:     # triggers dispense and drop tip if tip full or final transformant
    #         transfers_completed = (transformant // cell_tip_capacity) * cell_tip_capacity 
    #         req_transfers = transformant - transfers_completed + 1          # transfers required for this tip step (+1 to account for 0 indexing)

    #         for transfer in range(req_transfers):
    #             left.dispense(cells_vol + air_gap, TC_plate[t_array[transfers_completed + transfer]])
            
    #         left.touch_tip()

    #         tiprack300_1_count = iterate_tip_count(tiprack300_1_count, left)
    
    # # Transfer 3ul of plasmid to TC_plate
    # for plasmid in range(len(t_array)):
    #     right.pick_up_tip(tiprack20_1[tip_counter(tiprack20_1_count)])

    #     right.aspirate(DNA_vol, source_dict[plasmid][well_dict[plasmid]])

    #     right.dispense(DNA_vol, TC_plate[t_array[plasmid]])

    #     right.mix(DNA_vol, 10, TC_plate[t_array[plasmid]])

    #     right.touch_tip()

    #     tiprack20_1_count = iterate_tip_count(tiprack20_1_count, right)

    # 30 min incubation at 4c
    thermocycler.close_lid()
    thermocycler.set_block_temperature(4, hold_time_minutes = 30, block_max_volume = 55)
    
    # 30 second heatshock at 42c (NB: 30 rather than 45 due to slow ramping)
    thermocycler.set_block_temperature(42, hold_time_seconds = 30, block_max_volume = 55)
    
    # 90 second incubation at 4c
    thermocycler.set_block_temperature(4, hold_time_seconds = 90, block_max_volume = 55)
    
    # Add 100ul SOC per well
    thermocycler.open_lid()
    thermocycler.set_block_temperature(37)

    for transformant in range(len(t_array)):
        step = max(transformant, 1) % SOC_tube_capacity

        left.pick_up_tip(tiprack300_1[multi2single_tip_counter(tiprack300_1_count)])

        if step == 0:
            SOC_tube_count += 1                                         # changes SOC tube when tube has been emptied

        left.aspirate(100, tube_rack[tube_counter(SOC_tube_count)])

        left.dispense(100, TC_plate[t_array[transformant]])             # NB: uses one pipette per transfer to avoid cross contamination

        left.mix(2, 50, TC_plate[t_array[transformant]])

        tiprack300_1_count = iterate_tip_count(tiprack300_1_count, left)

    # 60 min incubation at 37c
    thermocycler.close_lid()
    thermocycler.set_block_temperature(37, hold_time_minutes = 60, block_max_volume = 150)

    # Plating 
    thermocycler.open_lid()
    # thermocycler.set_block_temperature(37)

    # for transformant in range(len(t_array)):
    #     step = transformant % plate_capacity

    #     if step == 0:
    #         plate = str(transformant // plate_capacity + 1)
    #         comment = "Place plate " + plate + " onto robot at deck position 4"
    #         protocol.pause(comment)
    #         flash(5)
        
    #     left.pick_up_tip(tiprack300_1[multi2single_tip_counter(tiprack300_1_count)])
        
    #     left.aspirate(plating_vol, TC_plate[t_array[transformant]])

    #     left.blow_out(plate_2[plate_counter_T(step)].top())
    #     left.touch_tip()

    #     tiprack300_1_count = iterate_tip_count(tiprack300_1_count, left)

    # Finish protocol
    thermocycler.deactivate()
    protocol.set_rail_lights(False)