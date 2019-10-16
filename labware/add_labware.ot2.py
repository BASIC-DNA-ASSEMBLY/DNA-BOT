from opentrons import labware, instruments, modules, robot 

if '4ti0131_trough-12' not in labware.list():
    custom_plate = labware.create('4ti0131_trough-12', grid=(12, 1), spacing=(9, 0), diameter=8.25, depth=39.22, volume=21000)
    robot.comment('4ti0131_trough-12 added to labware database')
else:
    robot.comment('4ti0131_trough-12 already exists in the OT2 database and has not been updated. Refer to Opentrons API and Custom_labware.csv')


if '4ti0136_96_deep-well' not in labware.list():
    custom_plate = labware.create('4ti0136_96_deep-well', grid=(12, 8), spacing=(9, 9), diameter=8.2, depth=39.15, volume=2200)
    robot.comment('4ti0136_96_deep-well added to labware database')
else:
    robot.comment('4ti0136_96_deep-well already exists in the OT2 database and has not been updated. Refer to Opentrons API and Custom_labware.csv')


if 'Nunc_Omnitray' not in labware.list():
    custom_plate = labware.create('Nunc_Omnitray', grid=(12, 8), spacing=(9, 9), diameter=5.5, depth=0, volume=0)
    robot.comment('Nunc_Omnitray added to labware database')
else:
    robot.comment('Nunc_Omnitray already exists in the OT2 database and has not been updated. Refer to Opentrons API and Custom_labware.csv')


if '4ti-0960_FrameStar' not in labware.list():
    custom_plate = labware.create('4ti-0960_FrameStar', grid=(12, 8), spacing=(9, 9), diameter=5.5, depth=15, volume=200)
    robot.comment('4ti-0960_FrameStar added to labware database')
else:
    robot.comment('4ti-0960_FrameStar already exists in the OT2 database and has not been updated. Refer to Opentrons API and Custom_labware.csv')


if 'tube-rack_E1415-1500' not in labware.list():
    custom_plate = labware.create('tube-rack_E1415-1500', grid=(6, 4), spacing=(20, 20), diameter=12, depth=38, volume=1500)
    robot.comment('tube-rack_E1415-1500 added to labware database')
else:
    robot.comment('tube-rack_E1415-1500 already exists in the OT2 database and has not been updated. Refer to Opentrons API and Custom_labware.csv')


if 'Eppendorf_30133366_plate_96' not in labware.list():
    custom_plate = labware.create('Eppendorf_30133366_plate_96', grid=(12,8), spacing=(9,9), diameter=5.5, depth=19, volume=250)
    robot.comment('Eppendorf_30133366_plate_96 added to labware database')
else:
    robot.comment('Eppendorf_30133366_plate_96 already exists in the OT2 database and has not been updated. Refer to Opentrons API and Custom_labware.csv')


if 'aluminium-block_4ti-0960_FrameStar' not in labware.list():
    custom_plate = labware.create('aluminium-block_4ti-0960_FrameStar', grid=(12, 8), spacing=(9, 9), diameter=5.5, depth=15, volume=200)
    robot.comment('aluminium-block_4ti-0960_FrameStar added to labware database')
else:
    robot.comment('aluminium-block_4ti-0960_FrameStar already exists in the OT2 database and has not been updated. Refer to Opentrons API and Custom_labware.csv')


for c in robot.commands():
    print(c)