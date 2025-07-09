def generate_final_assembly_dict_list(constructs_list, clips_df):
    ''' Runs the generate_final_assembly_dict function producing a list of 
    dictionaries of the assembly source and destination. Subsets the list of 
    assemblies so that 1) the maximum number of tips used does no to 96 (max assemblies 
    per plate). A sub dictionary is generated for each chunk and then function  
    returns a list of the resulting sub assembly dicts    
    This method ensures that subsequent assembly scripts do not reuse empty clip 
    wells while also keeping the correct destination well locations '''

    CONSTRUCT_COUNT_TOTAL = len(constructs_list)
    # ASSEMBLY_PLATE_COUNT = CONSTRUCT_COUNT_TOTAL // MAX_ASSEMBLIES_PER_PLATE + 1           # plus one to include final partially full plate ################# old version

    # # Error ########################### NO ERROR AS NO MAX ASSEMBLIES
    # if clips_df['number'].sum() > MAX_CLIPS_TOTAL:
    #     raise ValueError(
    #         'Number of CLIP reactions exceeds {}. Reduce number of constructs in construct.csv.'.format(MAX_CLIPS_TOTAL))

    # final_assembly_dict = generate_final_assembly_dict(constructs_list, clips_df) ###### old version (returns dictionary rather than keys and values)

    final_assembly_dict_keys, final_assembly_dict_values = generate_final_assembly_dict(constructs_list, clips_df)

    construct_count  = 0                                                                         # number of assemblies for current assembly script (<96)
    master_mix_tips = len(list({len(sublist[0]) for sublist in final_assembly_dict_values}))    # returns the length list of unique construct lengths in build - this represents the number of tips needed for MM transfer
    tip_count       = master_mix_tips

    assembly_subset_lower = 0                                                                   # set upper and lower bounds for subset of assemblies for a given plate
    assembly_subset_upper = 0

    keys    = []
    values  = []
    assembly_dict_list = []

    for i, construct in enumerate(final_assembly_dict_values):
        is_new_assembly_tips    = bool((tip_count + len(construct[0])) // (TIPS_PER_BOX * MAX_FINAL_ASSEMBLY_TIPRACKS)) # new assembly when out of tips (tips required exceeds cap)
        is_new_assembly_wells   = construct_count == MAX_ASSEMBLIES_PER_PLATE                                            # new assembly when over assembly limit (exceeds number of assembly plate wells)
        is_new_assembly_end     = i+1 == len(final_assembly_dict_values)                                                # new assembly when final construct reached
        is_new_assembly = is_new_assembly_tips or is_new_assembly_wells or is_new_assembly_end

        if is_new_assembly:
            assembly_subset_lower = assembly_subset_upper
            assembly_subset_upper = i

            if is_new_assembly_end:
                keys.append(tip_counter(construct_count))
                assembly_subset_upper = CONSTRUCT_COUNT_TOTAL            #

            values = final_assembly_dict_values[assembly_subset_lower:assembly_subset_upper]

            sub_assembly_dict = {keys[i]: values[i] for i in range(len(keys))}
            assembly_dict_list.append(sub_assembly_dict)                # generate and append sub_clip_dict to list - allows for multiple clip reactions

            keys    = []
            values  = []

            tip_count = master_mix_tips                                 # reset tips for new assembly and include tips for distribution of master mix
            construct_count = 0

        for tip in range(len(construct[0])):
            tip_count += 1
            print(i, '\ttip =', tip_count + tip)


        keys.append(tip_counter(construct_count))
        construct_count += 1 

    # assembly_dict_list_old = [] ############# old version - iterates by assembly and not by tip, subsets only by max assemblies not by max tips

    # for plate in range(ASSEMBLY_PLATE_COUNT):
    #     subset_lower = (plate * MAX_ASSEMBLIES_PER_PLATE)               # set upper and lower bounds for subset of assemblies for a given plate
    #     subset_upper = subset_lower + MAX_ASSEMBLIES_PER_PLATE

    #     if subset_upper > CONSTRUCT_COUNT_TOTAL:                               # set total number number of assemblies as upper bound if plate incomplete
    #         subset_upper = CONSTRUCT_COUNT_TOTAL

    #     # sub_assembly_df = constructs_list[subset_lower:subset_upper]
    #     # sub_assembly_dict = generate_final_assembly_dict(sub_assembly_df, clips_df)

    #     keys = [tip_counter(i) for i in list(range(subset_upper - subset_lower))]
    #     values = final_assembly_dict_values[subset_lower:subset_upper]
    #     sub_assembly_dict = {keys[i]: values[i] for i in range(len(keys))}
    #     # print(values, '\n')

    #     assembly_dict_list_old.append(sub_assembly_dict)                    # generate and append sub_clip_dict to list - allows for multiple clip reactions
    return assembly_dict_list