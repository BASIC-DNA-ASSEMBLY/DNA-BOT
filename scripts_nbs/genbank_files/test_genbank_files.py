from Bio.Seq import UnknownSeq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, SeqUtils
import genbank_files

# PARTS = genbank_files.generate_seqrecords()
RFP_SEQRECORD = SeqIO.read("p004_rfp.gb", "genbank")
RFP_SEQRECORD.upper()
IP_STR = "TCTGGTGGGTCTCTGTCC"
IS_STR = "GGCTCGGGAGACCTATCG"


def test_classes():
    unk_dna = UnknownSeq(20, alphabet=IUPAC.ambiguous_dna)
    seqrecord1 = SeqRecord(unk_dna, id="An unknown seq")
    test_assembly = genbank_files.SeqRecAssembly(seqrecord1)
    test_promoter = genbank_files.PromoterPart(seqrecord1.seq, "000")
    print(test_promoter.id)
    test_orf = genbank_files.OrfPart(seqrecord1.seq, "fp")
    print(test_orf.id)
    print(type(test_promoter))
    print(test_promoter)


def test_parts():
    for part in PARTS:
        print("\n", part, "\nLength of SeqRecord:", len(part.seq))
        

def debug_generate_parts():
    basic_rfp_seqrecord = next(
        (part for part in PARTS if part.id == "BASIC_RFP_ORF_v1.0")
        , None)
    print(basic_rfp_seqrecord.id)
    pos = 0
    while str(RFP_SEQRECORD.seq)[pos] == str(basic_rfp_seqrecord.seq)[pos]:
        pos += 1
    print(f"position of difference is {pos}")
    print(str(RFP_SEQRECORD.seq[:pos]))
    print(
        f"the length of the original is {len(RFP_SEQRECORD.seq)} while the modified version is {len(basic_rfp_seqrecord.seq)}")
    print(str(RFP_SEQRECORD.seq[pos:]))
    print(str(RFP_SEQRECORD.seq))
    print(str(basic_rfp_seqrecord.seq[pos:]))
    return RFP_SEQRECORD.seq == basic_rfp_seqrecord.seq, str(basic_rfp_seqrecord.seq)


def debug_iseq_loc():
    search = SeqUtils.nt_search(str(RFP_SEQRECORD.seq), IP_STR)
    return search