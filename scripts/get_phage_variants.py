# Import modules.
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


##############################################
#################### MAIN ####################
##############################################


def main():

    # Define I/O paths.
    concat_results = snakemake.input[0]
    phage_variants = snakemake.output[0]

    # Define header of data frame.
    results_header = ['query', 'subject', 'pid', 'qcov', 'scov']

    # Load dataframe and add header.
    results_df = pd.read_csv(concat_results, sep='\t')
    results_df.columns = results_header
    phages = list(set(results_df['query'].to_list()))

    # Parameters to define phage variants.
    min_pid = snakemake.params.min_pid[0]
    min_qcov = snakemake.params.min_qcov[0]
    min_scov = snakemake.params.min_scov[0]

    # Filter for phage variants.
    filt_variants_pid = (results_df['pid'] >= min_pid)
    filt_variants_qcov = (results_df['qcov'] >= min_qcov)
    filt_variants_scov = (results_df['scov'] >= min_scov)
    filt_variants = filt_variants_pid & filt_variants_qcov & filt_variants_scov
    results_df = results_df.loc[filt_variants]

    # Prepare data frame for transforming to pivot table.
    results_df = results_df[['query', 'subject']]
    results_df['value'] = 1

    # Transform data frame to pivot table.
    results_pivot = pd.pivot_table(results_df, \
    values='value', index=['query'], columns=['subject'], fill_value=0)

    # Column by columns get rows (phage names) that have value 1.
    variants = [getPhageVariants(results_pivot, column) for column in phages]
    # From phage names leave only unique - list of phage names is one variant.
    uqvariants = list(set(variants))

    # Save file as phage-variants.csv
    filerows = ['PV,phage\n']
    for i, phages in enumerate(uqvariants):
        nvariant = f'PV{i+1}'
        for phage in phages:
            row = f'{nvariant},{phage}\n'
            filerows.append(row)

    with open(phage_variants, 'w+') as f:
        f.write(''.join(filerows))


# Define functions.
def getPhageVariants(pivot_table, column):
    filt_variants = (pivot_table[column] == 1)
    phages = pivot_table.loc[filt_variants].index.to_list()
    phages = sorted(phages)
    return tuple(phages)


if __name__ == '__main__':
    main()
