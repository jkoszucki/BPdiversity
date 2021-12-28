# Import modules.
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from Bio import SeqIO
import multiprocessing as mp
import itertools

##############################################
#################### MAIN ####################
##############################################

##############################################
# Here instead of blast results blast #######
# results are given. #########################
##############################################

def main():
    # Define I/O paths.
    # NOTE: blast OR blast results.
    blast_results = snakemake.input[0]
    table_output = snakemake.output[0]
    phage_names = snakemake.params.phage_names
    extension = snakemake.params.extension
    phages_dir = snakemake.params.phages_dir

    # Remove lines starting with hash if results are in blast format.
    clearBLAST(blast_results)

    # Define data frames headers.
    blast_header = ['query', 'subject', 'pid','alilength','missmatch',\
                     'gapopen','qstart','qend','sstart','send','evalue',\
                     'bitscore']
    results_header = ['query', 'subject', 'pid', 'qcov', 'scov']

    # Load results of blast search and name columns.
    blast_df = pd.read_csv(blast_results, sep='\t', header=None)
    blast_df.columns = blast_header

    # Get fullnames of phages.
    blast_df['query'] = blast_df.apply(get_fullname, args=(['query', phage_names]), axis=1)
    blast_df['subject'] = blast_df.apply(get_fullname, args=(['subject', phage_names]), axis=1)


    print('Improve prefiltering of results! Concatenate results within the same regions!!!')
    # Prefilter megablast results.
    # Filter blast results - only significant hits.
    filt_pid = (blast_df['pid'] >= 0.75)
    filt_eval = (blast_df['evalue'] <= 10**-3)
    filt_bitscore = (blast_df['bitscore'] >= 100)
    filt_alignment = (blast_df['alilength'] >= 300)
    filt = filt_pid & filt_eval & filt_bitscore & filt_alignment
    blast_df = blast_df.loc[filt]

    # Calculate and add to data frame query and subject coverages.
    blast_df['qcov'] = blast_df.apply(getCoverage, args=(['qcov', extension, phages_dir]), axis=1)
    blast_df['scov'] = blast_df.apply(getCoverage, args=(['scov', extension, phages_dir]), axis=1)

    # Get list of unique phage names.
    phages = list(set(blast_df['query'].to_list()))

    # Prepare input for multiprocessing.
    processes_pool = mp.Pool()
    blast_df_list = [blast_df]*(len(phages)*len(phages))
    q1, q2 = zip(*list(itertools.product(phages, repeat=2)))
    df_q1_q2 = list(zip(blast_df_list, q1, q2))

    # Run analysis via multiproccessing.
    with mp.Pool() as pool:
        results_rows = pool.map(concatenateLocalAlignmentsUnpack, df_q1_q2)
    results_rows = [row[:-1] for row in results_rows]

    # blast phages results: all vs all.
    results_df = pd.DataFrame(results_rows)
    results_df.columns = results_header

    # Remove rows with qcov & scov equal 0.
    filt_empty_query = (results_df['qcov'] == 0.0)
    filt_empty_subject = (results_df['scov'] == 0.0)
    filt_empty_hits = filt_empty_query & filt_empty_subject

    results_df = results_df.loc[~filt_empty_hits]

    results_df.to_csv(table_output, sep='\t', index=False)


# Define functions.
def get_fullname(phage_row, column, phage_names):
    """Substitute phage ordinal number by phage full name."""
    for phage_name in phage_names:
        if phage_row[column] == int(phage_name.split('_')[0]):
            return phage_name


def getCoverage(phage_row, column, extension, phages_dir, qcov=True):
    """ Calculates coverage of queries/subjects.
        Returns numpy array of queries/subjects coverages from whole dataframe. """
    if column == 'qcov':
        end, start = 'qend', 'qstart'
        query = phage_row['query']
    elif column == 'scov':
        end, start = 'send', 'sstart'
        query = phage_row['subject']
    else:
        print('Error in getCoverage function')
        exit()

    phage_length = getLength(query, extension, phages_dir)
    seq_length = abs((phage_row[end] - phage_row[start]))

    seq_coverage = seq_length/phage_length
    seq_coverage = round(seq_coverage, 2)
    return seq_coverage


def getLength(query, extension, phages_dir):
    """ Calculates length of query/subject. """
    phage_path = Path(phages_dir, query + f'.{extension}')
    phage_record = list(SeqIO.parse(phage_path, extension))[0]
    return len(phage_record.seq)


def concatenateLocalAlignmentsUnpack(args):
    return concatenateLocalAlignments(*args)


def concatenateLocalAlignments(df, complete1, complete2, round_float=9):
    filt_query = df['query'].str.contains(complete1)
    filt_subject = df['subject'].str.contains(complete2)
    filt = filt_query & filt_subject

    temp_df = df.loc[filt]
    if not temp_df.empty:
        pid_mean = np.round(np.mean(np.array(temp_df['pid'].to_list())), round_float)
        qcov_sum = np.round(np.sum(np.array(temp_df['qcov'].to_list())), round_float)
        scov_sum = np.round(np.sum(np.array(temp_df['scov'].to_list())), round_float)
        qhits_locations = list(zip(temp_df['qstart'].to_list(), temp_df['qend'].to_list()))


    else:
        pid_mean = np.round(0, round_float)
        qcov_sum = np.round(0, round_float)
        scov_sum = np.round(0, round_float)
        qhits_locations = [(0,0)]

    return complete1, complete2, pid_mean, abs(qcov_sum), abs(scov_sum), qhits_locations


def clearBLAST(results):
    with open(results, 'r') as r:
        clearBLASTlines = [line for line in r.readlines() if line[0] != '#']
        clearBLAST = ''.join(clearBLASTlines)

    with open(results, 'w+') as  r:
        r.write(clearBLAST)


if __name__ == '__main__':
    main()
