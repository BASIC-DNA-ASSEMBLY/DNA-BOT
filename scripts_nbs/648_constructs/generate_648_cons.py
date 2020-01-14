import itertools


def main():
    constructs = make_cons()
    splitted_cons = constructs.split_cons(88)
    return print([len(build) for build in splitted_cons])


def make_cons():
    """Function for making the required 648 constructs

    """

    promoter_abrevs = ["105", "106", "101", "104"]
    promoters = [write_promoter(abrev) for abrev in promoter_abrevs]
    orf_abrevs = ["sfGFP", "RFP", "BFP"]
    orfs = [write_orf(abrev) for abrev in orf_abrevs]
    orf_perms = list(itertools.permutations(orfs))
    rbs_perms = list(itertools.product(range(1, 4), repeat=3))
    constructs = []
    for promoter in promoters:
        for orf_perm in orf_perms:
            for rbs_perm in rbs_perms:
                constructs.append(BasicCon("LMS", "BASIC_SEVA_37_CmR-p15A_v1.0", "LMP", promoter,
                                           assign_rbs_linker_utr(
                                               orf_perm[0]) + str(rbs_perm[0]), orf_perm[0],
                                           assign_rbs_linker_utr(
                                               orf_perm[1]) + str(rbs_perm[1]), orf_perm[1],
                                           assign_rbs_linker_utr(orf_perm[2]) + str(rbs_perm[2]), orf_perm[2]))
    return ConsCollection(*constructs)


def write_promoter(anderson_abrev):
    return "BASIC_L3S2P21_J23" + anderson_abrev + "_RiboJ"


def write_orf(orf_abrev):
    return "BASIC_" + orf_abrev + "_ORF_v1.0"


def assign_rbs_linker_utr(orf, utr_dict={write_orf("sfGFP"): "1", write_orf("BFP"): "2", write_orf("RFP"): "3"}):
    return "UTR" + utr_dict[orf] + "-RBS"


class BasicCon:
    """ Basic construct class containing part linker descriptions with each as a string.

    """

    def __init__(self, *parts_linkers):
        if not all((isinstance(part_linker, str) for part_linker in parts_linkers)):
            raise TypeError("A part or linker is not of type string")
        self.parts_linkers = list(parts_linkers)

    def __str__(self):
        return str(self.parts_linkers)


class ConsCollection:
    """ A class for organising and manipulating BasicCon classes

    """

    def __init__(self, *basic_cons):
        if not all((isinstance(basic_con, BasicCon) for basic_con in basic_cons)):
            raise TypeError("A basic construct is not of type basic construct")
        self.basic_cons = list(basic_cons)

    def split_cons(self, number):
        """Returns a nested list where self.basic_cons has been split into multiples of number or fewer
        e.g. number = 88, returns a nested list where each sub-list < 88 BasicCon objects.

        """
        cons_splitted = []
        i = 0
        while i < len(self.basic_cons):
            if len(self.basic_cons) - i > number:
                cons_splitted.append(self.basic_cons[i:i + number])
                i += number
            else:
                cons_splitted.append(self.basic_cons[i:])
                i = len(self.basic_cons)
        return cons_splitted


if __name__ == "__main__":
    main()
