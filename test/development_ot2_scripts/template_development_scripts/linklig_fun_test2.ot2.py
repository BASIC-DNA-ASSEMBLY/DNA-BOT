from opentrons import labware, instruments, modules, robot
# Robot class added.

linker_ligation_dict={"prefix_linker_wells": ["C4", "A8", "C8", "D8", "A4", "C4", "A8", "C8", "D8", "A4", "C4", "A8", "C8", "D8", "A4", "C4", "A8", "C8", "D8", "A4", "B4", "C4", "A8", "C8", "D8", "A4", "B4", "C4", "A8", "C8", "D8", "A4", "B4"], "suffix_linker_wells": ["E8", "G8", "H8", "G6", "F4", "E8", "G8", "H8", "G6", "F4", "E8", "G8", "H8", "G6", "F4", "E8", "G8", "H8", "G6", "E4", "F4", "E8", "G8", "H8", "G6", "E4", "F4", "E8", "G8", "H8", "G6", "E4", "F4"], "part_wells": ["A6", "B6", "D6", "F6", "A5", "A6", "B6", "D6", "F6", "A5", "A6", "B6", "D6", "F6", "A5", "A6", "B6", "D6", "F6", "A5", "B5", "A6", "B6", "D6", "F6", "A5",
    "B5", "A6", "B6", "D6", "F6", "A5", "B5"], "destination_wells": ["A1", "B1", "C1", "D1", "E1", "A2", "B2", "C2", "D2", "E2", "A3", "B3", "C3", "D3", "E3", "A4", "B4", "C4", "D4", "E4", "F4", "A5", "B5", "C5", "D5", "E5", "F5", "A6", "B6", "C6", "D6", "E6", "F6"], "part_vols": [1.3, 1.2, 1.1, 1.1, 1.6, 1.3, 1.2, 1.1, 1.1, 1.6, 1.3, 1.2, 1.1, 1.1, 1.6, 1.3, 1.2, 1.1, 1.1, 1.6, 1.2, 1.3, 1.2, 1.1, 1.1, 1.6, 1.2, 1.3, 1.2, 1.1, 1.1, 1.6, 1.2], "water_vols": [6.7, 6.8, 6.9, 6.9, 6.4, 6.7, 6.8, 6.9, 6.9, 6.4, 6.7, 6.8, 6.9, 6.9, 6.4, 6.7, 6.8, 6.9, 6.9, 6.4, 6.8, 6.7, 6.8, 6.9, 6.9, 6.4, 6.8, 6.7, 6.8, 6.9, 6.9, 6.4, 6.8]}


def linker_ligation(
    prefix_linker_wells,
	suffix_linker_wells,
	part_wells,
	destination_wells,
	part_vols,
	water_vols:

	# Constants
	CANDIDATE_TIPRACK_SLOTS=['3', '6', '9']
	TIPRACK_TYPE='tiprack-10ul'
	PIPETTE_TYPE='P10_Single'
	PIPETTE_MOUNT='right'
	SOURCE_PLATE_TYPE='4titude_96_skirt_PCR'
	SOURCE_PLATE_POSITION='2'
	DESTINATION_PLATE_TYPE='4ti-0960_FrameStar'
	DESTINATION_PLATE_POSITION='1'
	TUBE_RACK_TYPE='tube-rack_E1415-1500'
	TUBE_RACK_POSITION='5'
	MASTER_MIX_WELL='A1'
	WATER_WELL='A2'
	MASTER_MIX_VOLUME=20
	INITIAL_TIP='A1'

	# Tiprack slots
	total_tips=4 * len(part_wells)
	letter_dict={'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7}
	tiprack_1_tips = (13 - int(INITIAL_TIP[1:])) * 8 - letter_dict[INITIAL_TIP[0]]

	if total_tips > tiprack_1_tips:
		tiprack_num = 1 + (total_tips - tiprack_1_tips) // 96 + \
					(1 if (total_tips - tiprack_1_tips) % 96 > 0 else 0)
	else:
		tiprack_num = 1
	slots = CANDIDATE_TIPRACK_SLOTS[:tiprack_num]

	# Define labware
	tipracks = [labware.load(TIPRACK_TYPE, slot) for slot in slots]
	if PIPETTE_TYPE != 'P10_Single':
		print('Define labware must be changed to use', PIPETTE_TYPE)
		exit()

	pipette = instruments.P10_Single(mount=PIPETTE_MOUNT, tip_racks=tipracks)
	pipette.start_at_tip(tipracks[0].well(INITIAL_TIP))
	source_plate = labware.load(SOURCE_PLATE_TYPE, SOURCE_PLATE_POSITION)
	destination_plate = labware.load(
						DESTINATION_PLATE_TYPE, DESTINATION_PLATE_POSITION)
	tube_rack = labware.load(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
	master_mix = tube_rack.wells(MASTER_MIX_WELL)
	water = tube_rack.wells(WATER_WELL)

	# Transfers
	pipette.pick_up_tip()
	pipette.transfer(MASTER_MIX_VOLUME, master_mix,
					destination_plate.wells(destination_wells), new_tip='never')
	pipette.drop_tip()
	pipette.transfer(water_vols, water,
					destination_plate.wells(destination_wells), new_tip='always', blow_out=True)
	pipette.transfer(1, source_plate.wells(prefix_linker_wells),
					destination_plate.wells(destination_wells), new_tip='always', blow_out=True)
	pipette.transfer(1, source_plate.wells(suffix_linker_wells),
					destination_plate.wells(destination_wells), new_tip='always', blow_out=True)
	pipette.transfer(part_vols, source_plate.wells(part_wells),
					destination_plate.wells(destination_wells), new_tip='always', mix_after=(4, 5))


if __name__ == '__main__':
    linker_ligation(**linker_ligation_dict)

# Robot comments
for c in robot.commands():
    print(c)
