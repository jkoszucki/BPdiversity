# Load modules.
from Bio import SeqIO
from functools import partial

# Dedine functions.
def loadGenbank(genbankfile):
    return list(SeqIO.parse(genbankfile, 'genbank'))[0]

def getRBPs(record, key_phrases):
    rbps_seq = []
    for f in record.features:
        if f.type == 'CDS':
            for key_phrase in key_phrases:
                if key_phrase in f.qualifiers['product']:
                    rbps_seq.append(f.qualifiers['translation'])
                    break
    return *rbps_seq

# Define paths and params
genbank_files = snakemake.input
rbps_out = snakemake.output[0]
key_phrases = snakemake.params.rbp_key_phrases

records = list(map(loadGenbank, genbank_files))
map_getRBPs = partial(getRBPs, key_phrases=key_phrases)
rbps_sequences = list(map(map_getRBPs, records))
