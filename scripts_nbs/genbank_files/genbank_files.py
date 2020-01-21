from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from pathlib import Path
import csv


def main():
    # Define all the starting part and linker sequences as SeqRecord objects
    basic_parts = generate_seqrecords()
    for part in basic_parts:
        SeqIO.write(part, f"{part.id}.gb", "genbank")
    parts_linkers = basic_parts + generate_linkers()

    path_to_csv = Path.cwd().parents[1] / "examples" / "construct_csvs" / \
        "storch_et_al_cons" / "storch_et_al_cons.csv"
    seqrec_assemblies = []
    with open(path_to_csv, "r", newline="") as cons_csv:
        csv_reader = csv.reader(cons_csv)
        for ind, row in enumerate(csv_reader):
            if ind > 0:
                seqrec_assemblies.append(
                    SeqRecAssembly(
                        *(identify_basic_part(part, parts_linkers) for part in row[1:])
                    )
                )
    seqrec_assemblies = (seqrec_assembly.assemble_seq_record() for seqrec_assembly in seqrec_assemblies)
    SeqIO.write(seqrec_assemblies, "dnabot_constructs.gb", "genbank")


def identify_basic_part(target, candidate_list):
    candidates = (candidate.id for candidate in candidate_list)
    return candidate_list[candidates.index(target)]


class SeqRecAssembly:
    def __init__(self, *parts_linkers):
        self.parts_linkers = parts_linkers

    def assemble_seq_record(self):
        """Method returns the assembled SeqRecord

        """
        assembled_seq_record = SeqRecord(Seq(str()))
        for part_linker in self.parts_linkers:
            assembled_seq_record += part_linker.basic_slice()
        return assembled_seq_record

    @property
    def parts_linkers(self):
        return self._parts_linkers

    @parts_linkers.setter
    def parts_linkers(self, values):
        if not all(
                issubclass(value, SeqRecord) for value in values):
            raise TypeError(
                "Not all *parts_linkers are a subclass of SeqRecord."
            )
        if not all(
                callable(getattr(value, "basic_slice", None)) for value in values):
            raise TypeError(
                "Not all *parts_linkers have method: basic_slice."
            )
        self._parts_linkers = values


class BasicLinker(SeqRecord):
    def __init__(self, seq, id, **kwargs):
        super().__init__(seq=seq, id=id, **kwargs)
        self._linker_feature()

    def basic_slice(self):
        return self

    def _linker_feature(self):
        self.features.append(
            SeqFeature(
                type="misc_feature",
                location=FeatureLocation(2, len(self.seq), strand=+1),
                qualifiers={
                    "function": ["BASIC DNA assembly linker"],
                    "standard_name": [f"{self.id}"]
                }
            )
        )


class BasicPart(SeqRecord):
    IP_STR = "TCTGGTGGGTCTCTGTCC"
    IS_STR = "GGCTCGGGAGACCTATCG"

    def __init__(self, seq, id, **kwargs):
        super().__init__(seq=seq, id=id, **kwargs)
        self._ip_loc = self._find_iseq(
            self.IP_STR, "iP sequence"
        )
        self._is_loc = self._find_iseq(
            self.IS_STR, "iS sequence"
        )

    def basic_slice(self):
        if self._ip_loc < self._is_loc:
            return self[
                self._ip_loc + len(self.IP_STR):self._is_loc]
        elif self._ip_loc > self._is_loc:
            return self[self._ip_loc + len(self.IP_STR):] + self[:self._is_loc]
        else:
            raise ValueError("incorrect sequence used.")

    def _find_iseq(self, iseq_str, iseq_id="integrated sequence"):
        search_out = SeqUtils.nt_search(
            str(self.seq), iseq_str
        )
        if len(search_out) < 2:
            raise ValueError(f"{self.id} lacks {iseq_id}")
        return search_out[1]


class _AbrevPart(SeqRecord):
    def __init__(self, seq, abrev, prefix, suffix, version, **kwargs):
        id_from_abrev = prefix + abrev + suffix + "." + str(version)
        super().__init__(seq=seq, id=id_from_abrev, **kwargs)


class AbrevPromoter(_AbrevPart):
    def __init__(self, seq, abrev, prefix=None, suffix=None, version=1, **kwargs):
        if not prefix:
            prefix = "BASIC_L3S2P21_J23"
        if not suffix:
            suffix = "_RiboJ"
        super().__init__(seq, abrev, prefix, suffix, version, **kwargs)


class AbrevOrf(_AbrevPart):
    def __init__(self, seq, abrev, prefix=None, suffix=None, version=1, **kwargs):
        if not prefix:
            prefix = "BASIC_"
        if not suffix:
            suffix = "_ORF"
        super().__init__(seq, abrev, prefix, suffix, version, **kwargs)


def generate_seqrecords():
    parts = []
    parts.append(AbrevPromoter(
        Seq("CTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCCGTGCCTACTCTGGAAAATCTTTTACGGCTAGCTCAGTCCTAGGTACTATGCTAGCAGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA", IUPAC.unambiguous_dna), "105"    ))
    parts.append(AbrevPromoter(
        Seq("CTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCCGTGCCTACTCTGGAAAATCTTTTACGGCTAGCTCAGTCCTAGGTATAGTGCTAGCAGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA", IUPAC.unambiguous_dna), "106"
    ))
    parts.append(AbrevPromoter(
        Seq("CTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCCGTGCCTACTCTGGAAAATCTTTTACAGCTAGCTCAGTCCTAGGTATTATGCTAGCAGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA", IUPAC.unambiguous_dna), "101"
    ))
    parts.append(AbrevPromoter(
        Seq("CTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCCGTGCCTACTCTGGAAAATCTTTGACAGCTAGCTCAGTCCTAGGTATTGTGCTAGCAGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA", IUPAC.unambiguous_dna), "104"
    ))
    parts.append(AbrevOrf(
        Seq("ATGCGTAAAGGCGAAGAACTGTTCACGGGCGTAGTTCCGATTCTGGTCGAGCTGGACGGCGATGTGAACGGTCATAAGTTTAGCGTTCGCGGTGAAGGTGAGGGCGACGCGACCAACGGCAAACTGACCCTGAAGTTCATCTGCACCACCGGTAAACTGCCGGTGCCTTGGCCGACCTTGGTGACGACGTTGACGTATGGCGTGCAGTGTTTTGCGCGTTATCCGGACCACATGAAACAACACGATTTCTTCAAATCTGCGATGCCGGAGGGTTACGTCCAGGAGCGTACCATTTCCTTCAAGGATGATGGCACTTACAAAACTCGCGCAGAGGTTAAGTTTGAAGGTGACACGCTGGTCAATCGTATCGAATTGAAGGGTATCGACTTTAAAGAGGATGGTAACATTCTGGGCCATAAACTGGAGTATAACTTCAACAGCCATAATGTTTACATTACGGCAGACAAGCAAAAGAACGGCATCAAGGCCAATTTCAAGATTCGCCACAATGTTGAGGACGGTAGCGTCCAACTGGCCGACCATTACCAGCAGAACACCCCAATTGGTGACGGTCCGGTTTTGCTGCCGGATAATCACTATCTGAGCACCCAAAGCGTGCTGAGCAAAGATCCGAACGAAAAACGTGATCACATGGTCCTGCTGGAATTTGTGACCGCTGCGGGCATCACCCACGGTATGGACGAGCTGTATAAGCGTCCGTAA", IUPAC.unambiguous_dna), "sfGFP"
    ))
    parts.append(AbrevOrf(
        Seq("ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCTTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAA", IUPAC.unambiguous_dna), "RFP"
    ))
    parts.append(AbrevOrf(
        Seq("ATGTCCGAGTTGATCAAAGAGAACATGCATATGAAATTATATATGGAAGGCACTGTAGATAATCATCATTTTAAATGTACGTCGGAAGGCGAAGGTAAACCATATGAAGGTACGCAGACGATGCGCATCAAGGTGGTGGAGGGCGGTCCGCTGCCATTCGCTTTCGATATTTTAGCCACGAGCTTCCTCTACGGTTCTAAAACTTTCATCAATCACACGCAGGGTATTCCGGACTTCTTTAAACAGTCGTTCCCGGAGGGTTTCACCTGGGAACGCGTTACCACGTATGAAGATGGTGGTGTGCTTACGGCAACGCAGGACACGAGCCTTCAGGATGGGTGTTTGATTTACAACGTGAAAATTCGTGGTGTGAACTTCACGTCTAACGGCCCGGTGATGCAGAAAAAAACACTGGGTTGGGAAGCCTTTACCGAAACCCTGTATCCGGCGGACGGTGGCCTGGAAGGCCGTAATGATATGGCCTTGAAATTAGTCGGCGGTTCACACCTGATCGCGAACGCGAAAACAACCTATCGTAGTAAAAAACCAGCCAAAAACCTGAAAATGCCGGGCGTCTACTACGTAGACTACCGTCTGGAGCGCATTAAAGAGGCGAATAATGAAACCTATGTCGAGCAGCACGAAGTTGCGGTTGCACGCTATTGCGATCTGCCCAGCAAACTGGGCCACAAGCTTAATGGTAGCTAA", IUPAC.unambiguous_dna), "BFP"
    ))
    for part in parts:
        if isinstance(part, AbrevPromoter):
            part.features.append(SeqFeature(
                type="regulatory",
                location=FeatureLocation(0, len(part.seq), strand=+1),
                qualifiers={"note": ["promoter"], "standard_name": [part.id]}
            ))
        if isinstance(part, AbrevOrf):
            part.features.append(SeqFeature(
                type="CDS",
                location=FeatureLocation(0, len(part.seq), strand=+1),
                qualifiers={
                    "note": ["fluorescent reporter protein"],
                    "gene": [part.id],
                    "translation": str(part.translate().seq[:-1])
                }
            ))

    rfp_seqrecord = SeqIO.read(Path().cwd() / "p004_rfp.gb", "genbank")
    rfp_seqrecord.upper()
    mcherry_search = SeqUtils.nt_search(
        str(rfp_seqrecord.seq),
        "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCTTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAA"
    )
    parts_in_vectors = []
    for ind, part in enumerate(parts):
        seqrec = rfp_seqrecord[:mcherry_search[1]] + part + rfp_seqrecord[
            mcherry_search[1] + len(mcherry_search[0]):]
        seqrec.upper()
        seqrec.seq.alphabet = IUPAC.unambiguous_dna
        parts_in_vectors.append(
            BasicPart(
                seq=seqrec.seq,
                id=part.id,
                description=f"{part.id} BASIC part stored in AmpR pUC vector",
                features=seqrec.features
            )
        )
    backbone = SeqIO.read("benchling_BASIC_SEVA_37_CmR-p15A.1.gb", "genbank")
    backbone.upper()
    backbone.seq.alphabet = IUPAC.unambiguous_dna
    parts_in_vectors.append(BasicPart(
        seq=backbone.seq,
        id="BASIC_SEVA_37_CmR-p15A.1",
        description="backbone vector for BASIC DNA assembly containing \
        chloramphenicol resistance marker, p15A origin and \
            mScarlet counter selection marker",
        features=backbone.features
    )
    )
    for ind, part in enumerate(parts_in_vectors):
        parts_in_vectors[ind].name = part.id[:len(part.id)-2]
        parts_in_vectors[ind].annotations = {
            "organism": "Escherichia coli",
            "date": "17-JAN-2020",
            "accessions": [],
            "sequence_version": 1,
            "topology": "circular"
        }
    return parts_in_vectors


def generate_linkers():
    def generate_basic_linker(id, word_str, alphabet=IUPAC.unambiguous_dna):
        seq_str = "GG" + word_str.upper()
        basic_linker = BasicLinker(
            seq=Seq(seq_str, alphabet=alphabet),
            id=id 
        )
        return basic_linker
    
    linkers = []
    linkers.append(generate_basic_linker("UTR1-RBS1", "ctcgttgaacaccgtcTCAGGTAAGTATCAGTTGTAAatcacacaggactagtcc"))
    linkers.append(generate_basic_linker("UTR1-RBS2", "ctcgttgaacaccgtcTCAGGTAAGTATCAGTTGTAAaaagaggggaaatagtcc"))
    linkers.append(generate_basic_linker("UTR1-RBS3", "ctcgttgaacaccgtcTCAGGTAAGTATCAGTTGTAAaaagaggagaaatagtcc"))
    linkers.append(generate_basic_linker("UTR2-RBS1", "ctcgtgttactattggCTGAGATAAGGGTAGCAGAAAatcacacaggactagtcc"))
    linkers.append(generate_basic_linker("UTR2-RBS3", "ctcgtgttactattggCTGAGATAAGGGTAGCAGAAAaaagaggagaaatagtcc"))
    linkers.append(generate_basic_linker("UTR3-RBS1", "ctcggtatctcgtggtCTGACGGTAAAATCTATTGTAatcacacaggactagtcc"))
    linkers.append(generate_basic_linker("UTR3-RBS3", "ctcggtatctcgtggtCTGACGGTAAAATCTATTGTAaaagaggagaaatagtcc"))
    linkers.append(generate_basic_linker("LMP", "ctcgggtaagaactcgCACTTCGTGGAAACACTATTAtctggtgggtctctgtcc"))
    linkers.append(generate_basic_linker("LMS", "ctcgggagacctatcgGTAATAACAGTCCAATCTGGTGTaacttcggaatcgtcc"))
    return linkers


if __name__ == "__main__":
    main()
