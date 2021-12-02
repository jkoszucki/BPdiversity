# Import modules.
from pathlib import Path
import pandas as pd
import numpy as np


##############################################
#################### MAIN ####################
##############################################


def main():

    # Define I/O paths.
    concat_results = snakemake.input[0]
    output = snakemake.output[0]

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
    phages_similar_dict = {phage: getPhageVariants(results_pivot, phage) for phage in phages}
    # Get sorted phages that had the higher number of similar phages.
    phages_sorted_keys = sorted(phages_similar_dict, key=lambda k: len(phages_similar_dict[k]), reverse=True)

    ### Define variants ###
    variants_dict = phages_similar_dict.copy()
    # Iterate over each phage starting from the one that has highest number of similar phages.
    for phage in phages_sorted_keys:
        # If phage is still in dictionary go furher.
        if phage in variants_dict.keys():
            # Get similar phages to the given phage.
            phages_variant = list(phages_similar_dict[phage])
            # Update phage variants it the phage has already been asigned to the other variant.
            updated_phages_variant = phages_variant.copy()
            # Remove given phage from the list of similar phages.
            phages_variant.remove(phage)

            # Remove from dictionary of all phages ones that has already been assigned to the variant.
            for phage_already_assigned in phages_variant:
                # If phage that have been already assined cannot be removed form the dictionary
                # it means that it already has been assigned to the variant ....
                try:
                    del variants_dict[phage_already_assigned]
                except:
                # .... Thus it has to be removed from the phages of the new variant.
                    updated_phages_variant.remove(phage_already_assigned)

            variants_dict[phage] = tuple(updated_phages_variant)

        else: continue


    ### Switch representat phage (key) to a longest phage in the variant. ###
    variants_resorted_dict = {}
    # Sort phages accordingly to the length (the longest is first in list)
    for phage in variants_dict.keys():
        # Sorting list "in place". It will change the list in dictionary (the same list).
        phages_variant = variants_dict[phage]
        phages_variant = sorted(phages_variant, key=lambda phage: int(phage.split('_')[-1]) - int(phage.split('_')[-2]) + 1, reverse=True)
        # Get the longest phage.
        phage_longest = phages_variant[0]
        # Create new dictionary. Longest phage is the key and value (list) contains all of the variant phages.
        variants_resorted_dict[phage_longest] = phages_variant

    # Overwrite the variants dict with new dict, where the key is the longest phage.
    variants_dict = variants_resorted_dict

    # Prepare keys to sort from the "biggest" (most numerous) phage variant.
    phages_keys = sorted(variants_dict, key=lambda k: len(variants_dict[k]), reverse=True)

    # Prompt phage variants with representant (longest phage).
    for phage in phages_keys:
        print(f'Representant (longest phage): {phage}. Members of variant (sorted by lengths): ===>', ','.join(variants_dict[phage]))


    # Save file as phage-variants.csv
    filerows = ['PV,phage\n']
    for i, phage in enumerate(phages_keys):
        phages_variant = variants_dict[phage]
        nvariant = f'PV{i+1}'
        for phage_in_variant in phages_variant:
            row = f'{nvariant},{phage_in_variant}\n'
            filerows.append(row)

    with open(output, 'w+') as f:
        f.write(''.join(filerows))

    ### Extend phage variants table ###
    # Add information about representative phage of the variant.
    # Add information about capsule type of the phage.

    # Load metadata path.
    metadata = snakemake.params.metadata
    # Save file as phage-variants-extended.csv
    filerows = ['PV,representat,phage,K_locus,ST,K_locus_confidence\n']
    for i, phage in enumerate(phages_keys):
        phages_variant = variants_dict[phage]
        nvariant = f'PV{i+1}'
        for phage_in_variant in phages_variant:
            capsule_type, capsule_confidence = getCapsuleType(phage_in_variant, metadata)
            row = f'{nvariant},{phage},{phage_in_variant},{capsule_type},{ST},{capsule_confidence}\n'
            filerows.append(row)

    output_extended = Path(Path(output).parent, 'phage-variants-extended.csv')
    with open(output_extended, 'w+') as f:
        f.write(''.join(filerows))


# Define functions.
def getPhageVariants(pivot_table, column):
    filt_variants = (pivot_table[column] == 1)
    phages = pivot_table.loc[filt_variants].index.to_list()
    phages = sorted(phages)
    return tuple(phages)


def getCapsuleType(phage, metadata):
    medatada_df = pd.read_csv(metadata, sep=';', header='infer', na_values='-')
    genome_name = phage.split('.')[0]

    filt = (medatada_df['ID'] == genome_name)
    capsule_type = medatada_df.loc[filt]['K_locus'].values[0]
    capsule_confidence = medatada_df.loc[filt]['K_locus_confidence'].values[0]
    ST = medatada_df.loc[filt]['ST'].values[0]

    return str(capsule_type), str(capsule_confidence), str(ST)


if __name__ == '__main__':
    main()
