from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, SeqUtils
from Bio.SeqUtils.CheckSum import seguid
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from pathlib import Path
import csv
import datetime

DATE = datetime.datetime.now()
DEFAULT_ANNOTATIONS = {
                        "source": "synthetic construct",
                        "organism": "synthetic construct",
                        "taxonomy": ["other sequences", "artificial sequences"],
                        "date": DATE.strftime("%d-") + DATE.strftime("%b").upper() + DATE.strftime("-%Y"),
                        "accessions": [],
                        "sequence_version": 1,
                        "topology": "circular"
                    }


def main():
    # Define all the starting part and linker sequences as SeqRecord objects
    basic_parts = generate_parts()
    for part in basic_parts:
        SeqIO.write(part, f"{part.id}.gb", "genbank")
    parts_linkers = basic_parts + generate_linkers()

    path_to_csv = Path.cwd().parents[1] / "examples" / "construct_csvs" / \
        "storch_et_al_cons" / "storch_et_al_cons.csv"
    seqrec_assemblies = []
    with open(path_to_csv, "r", newline="") as cons_csv:
        csv_reader = csv.reader(cons_csv)
        for ind, row in enumerate(csv_reader):
            if ind > 0 and row[1]:
                row = [
                    part_linker_string for part_linker_string in row if part_linker_string]
                seqrec_assembly = SeqRecAssembly(
                    *(identify_basic_part(part, parts_linkers) for part in row[1:])
                )
                assembled_seqrec = seqrec_assembly.assemble_seq_record(
                    id=f"dnabot_{row[0]}",
                    name=f"dnabot_{row[0]}",
                    annotations=DEFAULT_ANNOTATIONS
                )
                assembled_seqrec.seq.alphabet = IUPAC.unambiguous_dna
                seqrec_assemblies.append(assembled_seqrec)
    SeqIO.write(seqrec_assemblies, "dnabot_constructs.gb", "genbank")
    clean_basic_genbank(
        "benchling_basic_seva_18_ampr-puc-1.gb", 
        "BASIC_SEVA_18_AmpR-pUC.1.gb",
        description="backbone vector for BASIC DNA assembly containing ampicillin resistance marker, pUC origin and mScarlet counter selection marker.",
        id="BASIC_SEVA_18_AmpR-pUC.1",
        name="BASIC_SEVA_18_AmpR-pUC"
    )


def identify_basic_part(target, candidate_list):
    candidates = [candidate.id for candidate in candidate_list]
    return candidate_list[candidates.index(target)]


class SeqRecAssembly:
    def __init__(self, *parts_linkers):
        self.parts_linkers = parts_linkers

    def assemble_seq_record(self, **kwargs):
        """Method returns the assembled SeqRecord. kwargs can be used for SeqRecord object data.

        """
        assembled_seq_record = SeqRecord(Seq(str()))
        for part_linker in self.parts_linkers:
            assembled_seq_record += part_linker.basic_slice()
        assembled_seq_record.id = seguid(assembled_seq_record.seq)
        assembled_seq_record.description = f"BASIC DNA Assembly of {[part_linker.id for part_linker in self.parts_linkers]}"
        if kwargs:
            for key, value in kwargs.items():
                setattr(assembled_seq_record, key, value)
        return assembled_seq_record

    @property
    def parts_linkers(self):
        return self._parts_linkers

    @parts_linkers.setter
    def parts_linkers(self, values):
        if not all(
                issubclass(type(value), SeqRecord) for value in values):
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
        self.kwargs = kwargs

    def basic_slice(self):
        returned_seqrec = SeqRecord(seq=self.seq, id=self.id, **self.kwargs)
        if self._ip_loc < self._is_loc:
            return returned_seqrec[
                self._ip_loc + len(self.IP_STR):self._is_loc]
        elif self._ip_loc > self._is_loc:
            return returned_seqrec[self._ip_loc + len(self.IP_STR):] + returned_seqrec[:self._is_loc]
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
        self.abrev = abrev
        super().__init__(seq, abrev, prefix, suffix, version, **kwargs)


def generate_parts():
    parts = []
    parts.append(AbrevPromoter(
        Seq("CTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCCGTGCCTACTCTGGAAAATCTTTTACGGCTAGCTCAGTCCTAGGTACTATGCTAGCAGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA", IUPAC.unambiguous_dna), "105"))
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
        Seq("ATGCGTAAAGGCGAAGAACTGTTCACGGGCGTAGTTCCGATTCTGGTCGAGCTGGACGGCGATGTGAACGGTCATAAGTTTAGCGTTCGCGGTGAAGGTGAGGGCGACGCGACCAACGGCAAACTGACCCTGAAGTTCATCTGCACCACCGGTAAACTGCCGGTGCCTTGGCCGACCTTGGTGACGACGTTGACGTATGGCGTGCAGTGTTTTGCGCGTTATCCGGACCACATGAAACAACACGATTTCTTCAAATCTGCGATGCCGGAGGGTTACGTCCAGGAGCGTACCATTTCCTTCAAGGATGATGGCACTTACAAAACTCGCGCAGAGGTTAAGTTTGAAGGTGACACGCTGGTCAATCGTATCGAATTGAAGGGTATCGACTTTAAAGAGGATGGTAACATTCTGGGCCATAAACTGGAGTATAACTTCAACAGCCATAATGTTTACATTACGGCAGACAAGCAAAAGAACGGCATCAAGGCCAATTTCAAGATTCGCCACAATGTTGAGGACGGTAGCGTCCAACTGGCCGACCATTACCAGCAGAACACCCCAATTGGTGACGGTCCGGTTTTGCTGCCGGATAATCACTATCTGAGCACCCAAAGCGTGCTGAGCAAAGATCCGAACGAAAAACGTGATCACATGGTCCTGCTGGAATTTGTGACCGCTGCGGGCATCACCCACGGTATGGACGAGCTGTATAAGCGTCCGTAA", IUPAC.unambiguous_dna), "sfGFP", suffix=""
    ))
    parts.append(AbrevOrf(
        Seq("ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCTTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAA", IUPAC.unambiguous_dna), "mCherry", suffix=""
    ))
    parts.append(AbrevOrf(
        Seq("ATGTCCGAGTTGATCAAAGAGAACATGCATATGAAATTATATATGGAAGGCACTGTAGATAATCATCATTTTAAATGTACGTCGGAAGGCGAAGGTAAACCATATGAAGGTACGCAGACGATGCGCATCAAGGTGGTGGAGGGCGGTCCGCTGCCATTCGCTTTCGATATTTTAGCCACGAGCTTCCTCTACGGTTCTAAAACTTTCATCAATCACACGCAGGGTATTCCGGACTTCTTTAAACAGTCGTTCCCGGAGGGTTTCACCTGGGAACGCGTTACCACGTATGAAGATGGTGGTGTGCTTACGGCAACGCAGGACACGAGCCTTCAGGATGGGTGTTTGATTTACAACGTGAAAATTCGTGGTGTGAACTTCACGTCTAACGGCCCGGTGATGCAGAAAAAAACACTGGGTTGGGAAGCCTTTACCGAAACCCTGTATCCGGCGGACGGTGGCCTGGAAGGCCGTAATGATATGGCCTTGAAATTAGTCGGCGGTTCACACCTGATCGCGAACGCGAAAACAACCTATCGTAGTAAAAAACCAGCCAAAAACCTGAAAATGCCGGGCGTCTACTACGTAGACTACCGTCTGGAGCGCATTAAAGAGGCGAATAATGAAACCTATGTCGAGCAGCACGAAGTTGCGGTTGCACGCTATTGCGATCTGCCCAGCAAACTGGGCCACAAGCTTAATGGTAGCTAA", IUPAC.unambiguous_dna), "mTagBFP2", suffix=""
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
                    "gene": [part.abrev],
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
        parts_in_vectors[ind].name = clean_seqrec_name(part)
        parts_in_vectors[ind].annotations = DEFAULT_ANNOTATIONS
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
    linkers.append(generate_basic_linker(
        "UTR1-RBS1", "ctcgttgaacaccgtcTCAGGTAAGTATCAGTTGTAAatcacacaggactagtcc"))
    linkers.append(generate_basic_linker(
        "UTR1-RBS2", "ctcgttgaacaccgtcTCAGGTAAGTATCAGTTGTAAaaagaggggaaatagtcc"))
    linkers.append(generate_basic_linker(
        "UTR1-RBS3", "ctcgttgaacaccgtcTCAGGTAAGTATCAGTTGTAAaaagaggagaaatagtcc"))
    linkers.append(generate_basic_linker(
        "UTR2-RBS1", "ctcgtgttactattggCTGAGATAAGGGTAGCAGAAAatcacacaggactagtcc"))
    linkers.append(generate_basic_linker(
        "UTR2-RBS3", "ctcgtgttactattggCTGAGATAAGGGTAGCAGAAAaaagaggagaaatagtcc"))
    linkers.append(generate_basic_linker(
        "UTR3-RBS1", "ctcggtatctcgtggtCTGACGGTAAAATCTATTGTAatcacacaggactagtcc"))
    linkers.append(generate_basic_linker(
        "UTR3-RBS3", "ctcggtatctcgtggtCTGACGGTAAAATCTATTGTAaaagaggagaaatagtcc"))
    linkers.append(generate_basic_linker(
        "LMP", "ctcgggtaagaactcgCACTTCGTGGAAACACTATTAtctggtgggtctctgtcc"))
    linkers.append(generate_basic_linker(
        "LMS", "ctcgggagacctatcgGTAATAACAGTCCAATCTGGTGTaacttcggaatcgtcc"))
    return linkers


def clean_basic_genbank(file_path, return_path, **kwargs):
    """Writes a cleaner genbank file for a BASIC part.

    Args:
        file_path: path to genbank file.
        return_path: path the cleaned genbank file is returned to.
        **kwargs: dict items used to update the seqrecord.

    """
    cleaned_genbank = SeqIO.read(file_path, "genbank")
    cleaned_genbank.name = clean_seqrec_name(cleaned_genbank)
    cleaned_genbank.annotations = DEFAULT_ANNOTATIONS
    if kwargs:
        for key, value in kwargs.items():
            setattr(cleaned_genbank, key, value)
    SeqIO.write(cleaned_genbank, return_path, "genbank")


def clean_seqrec_name(seqrec):
    return seqrec.id[:len(seqrec.id)-2]


if __name__ == "__main__":
    main()
