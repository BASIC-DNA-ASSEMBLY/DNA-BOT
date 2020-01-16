from Bio.Seq import UnknownSeq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import genbank_files

unk_dna = UnknownSeq(20, alphabet=IUPAC.ambiguous_dna)
seqrecord1 = SeqRecord(unk_dna, id="An unknown seq")
test_assembly = genbank_files.SeqRecAssembly(seqrecord1)
test_promoter = genbank_files.PromoterPart(seqrecord1.seq, "000")
print(test_promoter.id)
test_orf = genbank_files.OrfPart(seqrecord1.seq, "fp")
print(test_orf.id)
print(type(test_promoter))
print(test_promoter)
