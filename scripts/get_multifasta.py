from Bio import SeqIO

files = snakemake.input
multifasta = snakemake.output[0]
ext = snakemake.params.extension


def load_genbank(file):
    record = list(SeqIO.parse(file, str(ext)))[0]
    # record.description = ''
    # record.id = '_'.join(record.id.split('_')[:-2])
    return record

records = list(map(load_genbank, files))
n = SeqIO.write(records, multifasta, 'fasta')
