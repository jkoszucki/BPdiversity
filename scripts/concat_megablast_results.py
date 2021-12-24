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
# Here instead of mmseqs results blast #######
# results are given. #########################
##############################################

def main():
    # Define I/O paths.
    # NOTE: mmseqs OR blast results.
    mmseqs_results = snakemake.input[0]
    table_output = snakemake.output[0]

    # Remove lines starting with hash if results are in blast format.
    clearBLAST(mmseqs_results)

    # Define data frames headers.
    mmseqs_header = ['query', 'subject', 'pid','alilength','missmatch',\
                     'gapopen','qstart','qend','sstart','send','evalue',\
                     'bitscore']
    results_header = ['query', 'subject', 'pid', 'qcov', 'scov']

    # Load results of mmseqs search and name columns.
    mmseqs_df = pd.read_csv(mmseqs_results, sep='\t', header=None)
    mmseqs_df.columns = mmseqs_header

    # Make subject and query names pretty.
    # mmseqs_df['query'] = mmseqs_df.apply(lambda row: '_'.join(row['query'].split('_')[:4]), axis=1)
    # mmseqs_df['subject'] = mmseqs_df.apply(lambda row: '_'.join(row['subject'].split('_')[:4]), axis=1)


    print('Improve prefiltering of results! Concatenate results within the same regions!!!')
    # Prefilter megablast results.
    # Filter mmseqs results - only significant hits.
    filt_pid = (mmseqs_df['pid'] >= 0.75)
    filt_eval = (mmseqs_df['evalue'] <= 10**-3)
    filt_bitscore = (mmseqs_df['bitscore'] >= 100)
    filt = filt_pid & filt_eval & filt_bitscore
    mmseqs_df = mmseqs_df.loc[filt]

    # Calculate and add to data frame query and subject coverages.
    qcov, scov = getCoverage(mmseqs_df, qcov=True),  getCoverage(mmseqs_df, qcov=False)
    mmseqs_df['qcov'], mmseqs_df['scov'] = abs(qcov), abs(scov)

    # Get list of unique phage names.
    phages = list(set(mmseqs_df['query'].to_list()))

    # Prepare input for multiprocessing.
    processes_pool = mp.Pool()
    mmseqs_df_list = [mmseqs_df]*(len(phages)*len(phages))
    q1, q2 = zip(*list(itertools.product(phages, repeat=2)))
    df_q1_q2 = list(zip(mmseqs_df_list, q1, q2))

    # Run analysis via multiproccessing.
    with mp.Pool() as pool:
        results_rows = pool.map(concatenateLocalAlignmentsUnpack, df_q1_q2)
    results_rows = [row[:-1] for row in results_rows]

    # MMseqs phages results: all vs all.
    results_df = pd.DataFrame(results_rows)
    results_df.columns = results_header
    results_df.to_csv(table_output, sep='\t', index=False)


# Define functions.
def getCoverage(df, qcov=True):
    """ Calculates coverage of queries/subjects.
        Returns numpy array of queries/subjects coverages from whole dataframe. """
    if qcov:
        end, start = 'qend', 'qstart'
        queries = df['query']
    else:
        end, start = 'send', 'sstart'
        queries = df['subject']

    phage_lengths = [getLength(query) for query in queries]
    phage_lengths = np.array(phage_lengths)
    seq_lengths = np.array((df[end] - df[start]))
    seq_coverage = seq_lengths/phage_lengths
    return np.round(seq_coverage, 6)


def getLength(query):
    """ Calculates length of query/subject from its name (string). """
    length = str(query.split('_')[-4])
    # stop = int(query.split('_')[-1])
    return length


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
