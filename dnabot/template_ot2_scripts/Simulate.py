
#import sys; sys.path.append('c:/programdata/anaconda3/lib/site-packages')
#import sys; sys.path.append('c:/users/baldwin/appdata/roaming/python/python310/site-packages')
#from opentrons import simulate
#metadata = {'apiLevel': '2.19'}
##from clip_template_Thermocycler_Gen2_APIv2_19 import run
#run

import opentrons.simulate
protocol_file = open('dnabot\template_ot2_scripts\clip_template_Thermocycler_Gen2_APIv2_19.py')
opentrons.simulate.simulate(protocol_file)

#Labware
#plate = protocol.load_labware('corning_96_wellplate_360ul_flat', 1)
#tiprack_1 = protocol.load_labware('opentrons_96_tiprack_300ul', 2)
#pipettes
#p300 = protocol.load_instrument('p300_single_gen2', 'right', tip_racks=[tiprack_1])
#protocol.max_speeds['Z'] = 10
#complex commands
#p300.transfer(100, plate['A1'], plate['B1'], mix_before=(2,50), touch_tip=True, blow_out=True, blowout_location='destination well', new_tip='always') 
#block commands
#p300.pick_up_tip()
#p300.aspirate(100, plate['A2'])
#p300.dispense(100, plate['B2'])
#mix(repetitions, volume, location, rate)
#p300.mix(3, 50, plate['B2'], 0.5)
#p300.blow_out(plate['B2'].bottom(10))
#p300.return_tip()
#print out commands so we can check the simulation
#for line in protocol.commands(): 
#        print(line)
