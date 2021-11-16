from Bio import SeqIO

genbankfiles = snakemake.input
multifasta = snakemake.output[0]

def load_genbank(file):
    record = list(SeqIO.parse(file, 'genbank'))[0]
    record.description = ''
    record.id = '_'.join(record.id.split('_')[:-2])
    return record

records = list(map(load_genbank, genbankfiles))
n = SeqIO.write(records, multifasta, 'fasta')
