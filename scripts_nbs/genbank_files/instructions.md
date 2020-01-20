# Instructions for generating genbank files

- For each construct in storch_et_al_cons.csv generate a SeqRecAssembly object using part and linker SeqRecords.

## How are BASIC parts and linkers assembled?

- BASIC iP and iS sequences are identifed within BasicParts.
- The indexes of the intervening sequence are used to slice the SeqRecord.
  - If iS locus < iP locus, slice and add sequences from outside the intervening region between iS and iP. This ensures that the API considers all circular sequences. Print a warning that SeqFeatures may have been lost.
- Concatenate sliced BasicParts with linkers constant SeqRecords using the "+" operator.
- Return the resulting SeqRecord.
- Inherit from SeqRecord and add method/s.
- Constants for linkers.