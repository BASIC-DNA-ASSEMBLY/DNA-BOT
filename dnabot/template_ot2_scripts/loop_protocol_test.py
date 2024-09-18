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
    source_plate = protocol.load_labware('nest_96_wellplate_100ul_pcr_full_skirt', 1)
    destination_plate = protocol.load_labware('nest_96_wellplate_100ul_pcr_full_skirt', 2)
    tiprack_1 = protocol.load_labware('opentrons_96_tiprack_300ul', 4)
    
    #pipettes
    pipette_right = protocol.load_instrument('p300_single_gen2', 'right', tip_racks=[tiprack_1])

    pipette_right.well_bottom_clearance.aspirate = 1  # tip is 1 mm above well bottom
    pipette_right.well_bottom_clearance.dispense = 2 
    
    #setting pipette flow rates - acts as global parameter until redefined
    #sets the actual rate in ul/s - this needs to be appropriate for a specific pipette
    pipette_right.flow_rate.aspirate = 100 
    pipette_right.flow_rate.dispense = 100
    pipette_right.flow_rate.blow_out = 100
    
    for i in range(15):
        pipette_right.pick_up_tip()
        pipette_right.aspirate(100, source_plate.wells()[i])
        pipette_right.dispense(100, destination_plate.wells()[i])
        pipette_right.drop_tip()  


    for i in range(2):
        pipette_right.pick_up_tip()
        pipette_right.aspirate(100, source_plate.wells()[i])
        pipette_right.dispense(100, destination_plate.wells()[i])
        #slower aspiration rates will give less chance of cavitation
        #higher dispense rates gives better mixing 
        #slower final mix reduces chance of liquid remaining
        #rate=() sets rate as a fraction of the value defined in the global 'instrument_context.flow_rate' above
        pipette_right.aspirate(50, destination_plate.wells()[i].bottom(3), rate=0.5)
        pipette_right.dispense(50, destination_plate.wells()[i].bottom(1), rate=1)
        pipette_right.aspirate(50, destination_plate.wells()[i].bottom(2), rate=0.5)
        pipette_right.dispense(50, destination_plate.wells()[i].bottom(3), rate=0.75)
        #short delay allows liquid to settle in tip
        protocol.delay(seconds=0.5)
        #final slow aspirate and dispense to ensure good liquid recovery
        #push_out in final dispense adds specified volume beyond the stop
        pipette_right.aspirate(60, destination_plate.wells()[i].bottom(1), rate=0.25)
        pipette_right.dispense(60, destination_plate.wells()[i].bottom(3), rate=0.25, push_out=1)
        #moving to a specific height in the well before blowout
        pipette_right.move_to(destination_plate.wells()[i].top(-5), speed=5) # move to 2mm below the top of current well
        pipette_right.blow_out() #without argument blows out into current location
        # this version of touch tip sets the fraction of the tube radius for touching, so tip does not get bent
        # v_offset sets height negative number is below rim
        # speed = speed of tip movement so it does not crash around
        pipette_right.touch_tip(radius=0.9, v_offset=-5, speed=10) #without location blows out into current location
        pipette_right.drop_tip() #default is trash
   

    #print out commands so we can check the simulation
    #ONLY REQUIRED IN SIMULATION, NOT WHEN RUNNING ROBOT
    for line in protocol.commands(): 
        print(line)
