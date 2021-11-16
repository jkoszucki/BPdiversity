import pandas as pd

bacterial_metadata = snakemake.input.bacterial_metadata
phage_variants = snakemake.input.phage_variants[0]
output = snakemake.output[0]

#### ERROR PRONE SEPARATORS ###
bacterial_df = pd.read_csv(bacterial_metadata, sep=';')
bacterial_df = bacterial_df[['ID', 'ST', 'K_locus']]

phage_variants = pd.read_csv(phage_variants, sep=',')
phage_variants['ID'] = phage_variants['phage'].str.split('.').str[0]
phage_variants['value'] = 1
phage_variants.drop(columns=['phage'], inplace=True)

phage_variants_pivot = pd.pivot_table(phage_variants, \
values='value', index=['ID'], columns=['PV'], fill_value=0)

metadata_df = pd.merge(bacterial_df, phage_variants_pivot, on='ID', how='outer')
metadata_df = metadata_df.fillna(0)
metadata_df.to_csv(output, sep=',', index=False)
