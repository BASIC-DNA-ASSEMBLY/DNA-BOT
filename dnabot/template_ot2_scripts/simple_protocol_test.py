from opentrons import protocol_api
#may not be required, but can be used in code if specified 
# - needs to be done for any other packages used directly in code
import numpy as np

#metadata
metadata = {
     'apiLevel': '2.19',
     'protocolName': 'Simple code test',
     'description': 'use to evaluate simple pipetting and coding commands'
}

#Function that executes program
def run(protocol: protocol_api.ProtocolContext):
    #Labware
    #Use generic opentrons labware definitions if simulating
    #setting up custom labware locations for simulations is not straightforward
    #WHEN RUNNING IN THE LAB THESE MUST BE SUBSTITUTED FOR ACTUAL LABWARE BEING USED
    source_plate = protocol.load_labware('nest_96_wellplate_100ul_pcr_full_skirt', 1)
    destination_plate = protocol.load_labware('nest_96_wellplate_100ul_pcr_full_skirt', 2)
    tiprack_1 = protocol.load_labware('opentrons_96_tiprack_300ul', 4)
    
    #pipettes
    pipette_right = protocol.load_instrument('p300_single_gen2', 'right', tip_racks=[tiprack_1])
    
    #setting well clearances - acts as global parameter until redefined
    #when setting need to consider volume and well characteristics
    #important that liquid makes contact with well base, smaller volumes will require tighter clearance
    #if liquid is already present in the well, will want to pipette into the liquid
    #if tip touches base of well this will block aspirate and dispense functions
    pipette_right.well_bottom_clearance.aspirate = 1  # tip is 1 mm above well bottom
    pipette_right.well_bottom_clearance.dispense = 2 
    
    #setting pipette flow rates - acts as global parameter until redefined
    #sets the actual rate in ul/s - this needs to be appropriate for a specific pipette
    #Default flow rates - can increase these up to 2x or 3x (may have issues at higher speed)
    #P20 Gen2 7.56 ul/s
    #P300 Gen2 92.86 ul/s
    #P1000 Gen2 274.7 ul/s
    pipette_right.flow_rate.aspirate = 100 
    pipette_right.flow_rate.dispense = 100
    pipette_right.flow_rate.blow_out = 200
    
    #complex commands
    pipette_right.transfer(
        100, 
        source_plate['A1'].bottom(2), 
        destination_plate['B1'].bottom(5),
        mix_before=(2,50), 
        touch_tip=True, 
        blow_out=True,
        blowout_location="destination well", 
        new_tip='always',
        ) 
    
    # block commands
    # Set of commands below does a single transfer followed by a mix
    # Defining each step in turn enables much finer control of pipetting actions
    # VALUES SHOWN BELOW HAVE NOT BEEN VALIDATED IN TESTING - that's your job!
    pipette_right.pick_up_tip()
    pipette_right.aspirate(100, source_plate['A2'].bottom(1), rate=0.5)
    pipette_right.dispense(100, destination_plate['B2'].bottom(2), rate=0.25)
    # mix-block commands - can separately define rates and well heights for each step
    # Don't have to aspirate and dispense from same height
    # slower aspiration rates will give less chance of cavitation
    # higher dispense rates gives better mixing 
    # slower final mix reduces chance of liquid remaining
    # rate=() sets rate as a fraction of the value defined in the global 'instrument_context.flow_rate' above
    pipette_right.aspirate(50, destination_plate['B2'].bottom(3), rate=0.5)
    pipette_right.dispense(50, destination_plate['B2'].bottom(1), rate=1)
    pipette_right.aspirate(50, destination_plate['B2'].bottom(1), rate=0.5)
    pipette_right.dispense(50, destination_plate['B2'].bottom(3), rate=0.75)
    # short delay allows any excess liquid to settle in tip
    protocol.delay(seconds=0.5)
    # final slow aspirate and dispense to ensure good liquid recovery
    # push_out in final dispense adds specified volume beyond the stop
    pipette_right.aspirate(60, destination_plate['B2'].bottom(1), rate=0.25)
    pipette_right.dispense(60, destination_plate['B2'].bottom(3), rate=0.25, push_out=1)
    # SLOWLY moving to a specific height in the well before blowout
    pipette_right.move_to(destination_plate['B2'].top(-5), speed=5) # move to 2mm below the top of current well
    # OR position can be specified in the .blow_out command, cannot include a rate OR speed argument.
    # if no argument provided .blow_out() will blow out to current location
    pipette_right.blow_out(destination_plate['B2'].top(-5))
    # this version of touch tip sets the fraction of the tube radius for touching, so tip does not get bent
    # v_offset sets height negative number is below rim
    # speed = speed of tip movement so it does not crash around
    pipette_right.touch_tip(radius=0.9, v_offset=-5, speed=10)
    pipette_right.drop_tip()
   

    #print out commands so we can check the simulation
    #ONLY REQUIRED IN SIMULATION, NOT WHEN RUNNING ROBOT
    for line in protocol.commands(): 
        print(line)
