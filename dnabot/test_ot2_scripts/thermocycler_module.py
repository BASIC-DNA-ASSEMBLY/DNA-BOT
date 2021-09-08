#thermocycler module
from opentrons import  simulate
protocol = simulate.get_protocol_api('2.8')

#thermocycler block(destination_plate), default slot:7,8,10,11
tc_mod = protocol.load_module('Thermocycler Module')
destination_plate = tc_mod.load_labware('nest_96_wellplate_100ul_pcr_full_skirt')


 
#open lid and add reagent directly to thermocycler block
tc_mod.open_lid()

#transfer code here

 #run clip reaction in thermocycler block
tc_mod.close_lid()
tc_mod.set_lid_temperature(37)#set the lid temperature to 37˚C before running clip reaction
profile = [
             {'temperature': 37, 'hold_time_minutes': 2},
             {'temperature': 20, 'hold_time_minutes': 1}]#run digestion and ligation for 20 cycles
tc_mod.execute_profile(steps=profile, repetitions=20,block_max_volume=30)
tc_mod.set_block_temperature(60, hold_time_minutes=10, block_max_volume=30)#60˚C for 10 minutes
tc_mod.open_lid()#complete clip reaction
 
for c in protocol.commands():
   print(c)