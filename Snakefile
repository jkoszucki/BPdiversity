from pathlib import Path

configfile: 'config.yaml'

# Distances between bacterial isolates (mash table).
MASHTABLE = config['input'][0]
# Bacterial metadata (capsule & ST)
BACTERIAL_METADATA = config['input'][1]
# Directory to (pro)phages in genbank format.
PHAGES_DIR = config['input'][2]
# Output directory.
OUTPUT_DIR = config['output_dir'][0]

# Extension of prophage files in genbank format.
EXTENSION, = config.get('files_extension', 'fasta')
# Number of threads to be used.
THREADS = config.get('threads', 8)
# Lineage cutoff to be used (eg, 0.01 or 0.004)
CUTOFF = config.get('lineage_cutoff', 0)[0]
phages, = glob_wildcards(Path(PHAGES_DIR, "{genome}." + EXTENSION))

# Prompt given input.
print(f'Loaded bacterial distances (patristic) matrix: {MASHTABLE} with cutoff {CUTOFF}')
print(f'Running pipeline on genomes from folder: {PHAGES_DIR}')
print(f'With extension: {EXTENSION}')
print(f'Metadata loaded from: {BACTERIAL_METADATA}')
print(f'Rusults in: {OUTPUT_DIR}')

# Check input.
infiles = [MASHTABLE, phages, BACTERIAL_METADATA]
messages = ['mash table', 'phage genomes', 'bacterial metadata']

for infile, msg in zip(infiles, messages):
    if not infile:
        # Prompt missing input.
        print(f'No {msg} given. Aborting!')
        exit()


rule target:
    input:
        Path(OUTPUT_DIR, 'phages', 'phage-variants.csv'),
        # Path(OUTPUT_DIR, 'rbps', 'rbp-variants.csv'),
        # Path(OUTPUT_DIR, f'metadata.csv'),
        # Path(OUTPUT_DIR, f'phyloheatmap2-{CUTOFF}.pdf')


rule multifasta:
    input: expand(Path(PHAGES_DIR, "{phage}." + EXTENSION), phage=phages)
    output: Path(OUTPUT_DIR, "phages", "phages.mf")
    params: extension=EXTENSION,
    conda: 'scripts/getmultifasta.yaml'
    script: 'scripts/get_multifasta.py'


rule makeDB:
    input: rules.multifasta.output
    output: directory(Path(OUTPUT_DIR, "phages", "blastdb"))
    shell: 'makeblastdb -in {input} -dbtype nucl -out {output}/phagesdb'


rule megablast:
    input:
        query=rules.multifasta.output,
        db=rules.makeDB.output
    output: Path(OUTPUT_DIR, "phages", "blast-results.txt")
    shell:
        'blastn -task megablast -outfmt 7 -query {input.query} -db {input.db}/phagesdb >> {output} '


rule concatmegablast:
    input: rules.megablast.output
    output: Path(OUTPUT_DIR, "phages", "blast-results-concat.txt")
    # Prefiltering of megablast results! Hardcoded parameters at the moment!
    params:
        phage_names=phages,
        extension=EXTENSION,
        phages_dir=PHAGES_DIR

        # min_pid=config['prefilt']['min_pid'],
        # min_scov=config['prefilt']['min_subject_cov'],
        # min_qcov=config['prefilt']['min_query_cov']
    conda: 'scripts/phagevariants.yaml'
    script: 'scripts/concat_megablast_results.py'


rule phageVariants:
    input: rules.concatmegablast.output
    output: Path(OUTPUT_DIR, "phages", "phage-variants.csv")
    params:
        min_pid=config['phage_variants']['min_pid'],
        min_scov=config['phage_variants']['min_subject_cov'],
        min_qcov=config['phage_variants']['min_query_cov'],
        metadata=BACTERIAL_METADATA
    conda: 'scripts/phagevariants.yaml'
    script: 'scripts/get_phage_variants.py'


# rule rbpProteins:
#     input: expand(Path(PHAGES_DIR, "{phage}." + EXTENSION), phage=phages)
#     output: Path(OUTPUT_DIR, 'rbps', 'rbp-proteins.csv')
#     params:
#         rbp_key_words=config['receptor_binding_proteins']['rbp_key_phrases']
#     conda: 'scripts/phagevariants.yaml'
#     scripts: 'scripts/get_rbp_proteins.py'

# rule rbpVariants:
# aa_pid=config['receptor_binding_proteins']['aa_pid'],
# aa_cov=config['receptor_binding_proteins']['aa_cov']


rule results:
    input:
        bacterial_metadata=BACTERIAL_METADATA,
        phage_variants=rules.phageVariants.output
    output: Path(OUTPUT_DIR, f'metadata.csv')
    conda: 'scripts/phagevariants.yaml'
    script: 'scripts/get_results.py'


rule treeAndHeatmap:
    input:
        mashtable=MASHTABLE,
        results=rules.results.output
    params: CUTOFF
    output: Path(OUTPUT_DIR, f'phyloheatmap2-{CUTOFF}.pdf')
    conda: 'scripts/treeandheatmap.yaml'
    script: 'scripts/tree_and_heatmap.R'
