from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pathlib


def main():
    # Define all the starting part and linker sequences as SeqRecords

    """For each construct in storch_et_al_cons.csv generate a SeqRecAssembly object using part and linker SeqRecords. 
    Generate a dict which maps csv strings to SeqRecord objects. """

    """How are BASIC parts and linkers assembled?
    
    - BASIC iP and iS sequences are identifed within BasicParts. 
    - The indexes of the intervening sequence are used to slice the SeqRecord.
         - If iS locus < iP locus, slice and add sequences from outside the intervening region
         between iS and iP. This ensures that the API considers all circular sequences. Print
         a warning that SeqFeatures may have been lost.
    - Concatenate sliced BasicParts with linkers constant SeqRecords using the "+" operator.
    - Return the resulting SeqRecord.

    Inherit from SeqRecord and add method/s.
    constants for linkers
    
    """
    parts = generate_seqrecords()
    for part in parts:
        SeqIO.write(part, f"{part.id}.gb", "genbank")
    
    # seqrec_assemblies = (seqrec_assembly.assemble_seq_record() for seqrec_assembly in seqrec_assemblies)
    # SeqIO.write(seqrec_assemblies, "dnabot_constructs.gb", "genbank")


class SeqRecAssembly:
    def __init__(self, *parts_linkers):
        self.parts_linkers = parts_linkers

    def assemble_seq_record(self):
        """assemble function returns the 

        """
        pass


class _DnaPart(SeqRecord):
    def __init__(self, seq, abrev, prefix, suffix):
        basic_id = prefix + abrev + suffix
        super().__init__(seq=seq, id=basic_id)


class PromoterPart(_DnaPart):
    def __init__(self, seq, abrev, prefix=None, suffix=None):
        if not prefix:
            prefix = "BASIC_L3S2P21_J23"
        if not suffix:
            suffix = "_RiboJ.1"
        super().__init__(seq, abrev, prefix, suffix)


class OrfPart(_DnaPart):
    def __init__(self, seq, abrev, prefix=None, suffix=None):
        if not prefix:
            prefix = "BASIC_"
        if not suffix:
            suffix = "_ORF.1"
        super().__init__(seq, abrev, prefix, suffix)


def generate_seqrecords():
    parts = []
    parts.append(PromoterPart(
        Seq("CTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCCGTGCCTACTCTGGAAAATCTTTTACGGCTAGCTCAGTCCTAGGTACTATGCTAGCAGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA"
        , IUPAC.unambiguous_dna), "105"
    ))
    parts.append(PromoterPart(
        Seq("CTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCCGTGCCTACTCTGGAAAATCTTTTACGGCTAGCTCAGTCCTAGGTATAGTGCTAGCAGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA"
        , IUPAC.unambiguous_dna), "106"
    ))
    parts.append(PromoterPart(
        Seq("CTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCCGTGCCTACTCTGGAAAATCTTTTACAGCTAGCTCAGTCCTAGGTATTATGCTAGCAGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA"
        , IUPAC.unambiguous_dna), "101"
    ))
    parts.append(PromoterPart(
        Seq("CTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCCGTGCCTACTCTGGAAAATCTTTGACAGCTAGCTCAGTCCTAGGTATTGTGCTAGCAGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA"
        , IUPAC.unambiguous_dna), "104"
    ))
    parts.append(OrfPart(
        Seq("ATGCGTAAAGGCGAAGAACTGTTCACGGGCGTAGTTCCGATTCTGGTCGAGCTGGACGGCGATGTGAACGGTCATAAGTTTAGCGTTCGCGGTGAAGGTGAGGGCGACGCGACCAACGGCAAACTGACCCTGAAGTTCATCTGCACCACCGGTAAACTGCCGGTGCCTTGGCCGACCTTGGTGACGACGTTGACGTATGGCGTGCAGTGTTTTGCGCGTTATCCGGACCACATGAAACAACACGATTTCTTCAAATCTGCGATGCCGGAGGGTTACGTCCAGGAGCGTACCATTTCCTTCAAGGATGATGGCACTTACAAAACTCGCGCAGAGGTTAAGTTTGAAGGTGACACGCTGGTCAATCGTATCGAATTGAAGGGTATCGACTTTAAAGAGGATGGTAACATTCTGGGCCATAAACTGGAGTATAACTTCAACAGCCATAATGTTTACATTACGGCAGACAAGCAAAAGAACGGCATCAAGGCCAATTTCAAGATTCGCCACAATGTTGAGGACGGTAGCGTCCAACTGGCCGACCATTACCAGCAGAACACCCCAATTGGTGACGGTCCGGTTTTGCTGCCGGATAATCACTATCTGAGCACCCAAAGCGTGCTGAGCAAAGATCCGAACGAAAAACGTGATCACATGGTCCTGCTGGAATTTGTGACCGCTGCGGGCATCACCCACGGTATGGACGAGCTGTATAAGCGTCCGTAA"
        , IUPAC.unambiguous_dna), "sfGFP"
    ))
    parts.append(OrfPart(
        Seq("ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCTTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAA"
        , IUPAC.unambiguous_dna), "RFP"
    ))
    parts.append(OrfPart(
        Seq("ATGTCCGAGTTGATCAAAGAGAACATGCATATGAAATTATATATGGAAGGCACTGTAGATAATCATCATTTTAAATGTACGTCGGAAGGCGAAGGTAAACCATATGAAGGTACGCAGACGATGCGCATCAAGGTGGTGGAGGGCGGTCCGCTGCCATTCGCTTTCGATATTTTAGCCACGAGCTTCCTCTACGGTTCTAAAACTTTCATCAATCACACGCAGGGTATTCCGGACTTCTTTAAACAGTCGTTCCCGGAGGGTTTCACCTGGGAACGCGTTACCACGTATGAAGATGGTGGTGTGCTTACGGCAACGCAGGACACGAGCCTTCAGGATGGGTGTTTGATTTACAACGTGAAAATTCGTGGTGTGAACTTCACGTCTAACGGCCCGGTGATGCAGAAAAAAACACTGGGTTGGGAAGCCTTTACCGAAACCCTGTATCCGGCGGACGGTGGCCTGGAAGGCCGTAATGATATGGCCTTGAAATTAGTCGGCGGTTCACACCTGATCGCGAACGCGAAAACAACCTATCGTAGTAAAAAACCAGCCAAAAACCTGAAAATGCCGGGCGTCTACTACGTAGACTACCGTCTGGAGCGCATTAAAGAGGCGAATAATGAAACCTATGTCGAGCAGCACGAAGTTGCGGTTGCACGCTATTGCGATCTGCCCAGCAAACTGGGCCACAAGCTTAATGGTAGCTAA"
        , IUPAC.unambiguous_dna), "BFP"
    ))
    for part in parts:
        if isinstance(part, PromoterPart):
            part.features.append(SeqFeature(
                type="regulatory",
                location=FeatureLocation(0, len(part.seq), strand=+1),
                qualifiers={"note": ["promoter"], "standard_name": [part.id]}
                ))
        if isinstance(part, OrfPart):
            part.features.append(SeqFeature(
                type="CDS",
                location=FeatureLocation(0, len(part.seq), strand=+1),
                qualifiers={
                    "note": ["fluorescent reporter protein"],
                    "gene": [part.id],
                    "translation": str(part.translate().seq[:-1])
                    }
                ))
    rfp_seqrecord = SeqIO.read(pathlib.Path().cwd() / "p004_rfp.gb", "genbank")
    rfp_seqrecord.upper()
    mcherry_search = SeqUtils.nt_search(
        str(rfp_seqrecord.seq),
        "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCTTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAA"
        )
    parts_in_vectors = []
    for ind, part in enumerate(parts):
        parts_in_vectors.append(rfp_seqrecord[:mcherry_search[1]] + part + rfp_seqrecord[
        mcherry_search[1] + len(mcherry_search[0]):])
        parts_in_vectors[ind].id = part.id
        parts_in_vectors[ind].name = part.id[:len(part.id)-2]      
        parts_in_vectors[ind].description = f"{part.id} BASIC part stored in AmpR pUC vector"
        parts_in_vectors[ind].annotations = {
            "organism": "Escherichia coli",
            "date": "17-JAN-2020",
            "accessions": [],
            "sequence_version": 1,
            "topology": "circular"
        }
    return parts_in_vectors


if __name__ == "__main__":
    main()