# pipeline to generate read counts and phenotype scores directly from gzipped sequencing data

import os
import sys
import gzip
import multiprocessing
import fnmatch
import glob
import argparse


# CONSTANTS

ACCEPTED_FILE_TYPES = {'*.fastq.gz': 'fqgz',
                       '*.fastq': 'fq',
                       '*.fq': 'fq',
                       '*.fa': 'fa',
                       '*.fasta': 'fa',
                       '*.fna': 'fa'}

TEST_LINES = 10000


def parse_arguments(parser=None):
    if not parser:
        parser = argparse.ArgumentParser(description='Process raw sequencing data from screens to counts '
                                                     'files in parallel')
    parser.add_argument('Library_Fasta', help='Fasta file of expected library reads.')
    parser.add_argument('Out_File_Path', help='Directory where output files should be written.')
    parser.add_argument('Seq_File_Names', nargs='+',
                        help='Name(s) of sequencing file(s). Unix wildcards can be used to select multiple files '
                             'at once. The script will search for all *.fastq.gz, *.fastq, and *.fa(/fasta/fna) '
                             'files with the given wildcard name.')
    parser.add_argument('-p', '--processors', type=int, default=1)
    parser.add_argument('--trim_start', type=int)
    parser.add_argument('--trim_end', type=int)
    parser.add_argument('--test', action='store_true', default=False,
                        help='Run the entire script on only the first %d reads of each file. Be sure to delete or '
                             'move all test files before re-running script as they will not be overwritten.'
                             % TEST_LINES)
    arguments = parser.parse_args()

    return arguments


def parallel_seq_file_to_counts_parallel(fastq_gz_file_name_list, fasta_file_name_list, count_file_name_list,
                                         process_pool, library_fasta, start_index=None, stop_index=None, test=False):
    """Sequence File to Trimmed Fasta Functions"""
    if len(fastq_gz_file_name_list) != len(fasta_file_name_list):
        raise ValueError('In and out file lists must be the same length')

    arglist = list(
        zip(fastq_gz_file_name_list, fasta_file_name_list, count_file_name_list,
            [library_fasta] * len(fasta_file_name_list),
            [start_index] * len(fasta_file_name_list), [stop_index] * len(fasta_file_name_list),
            [test] * len(fasta_file_name_list)))

    reads_per_file = process_pool.map(seq_file_to_counts_wrapper, arglist)

    return list(zip(count_file_name_list, reads_per_file))


def seq_file_to_counts_wrapper(arg):
    return seq_file_to_counts(*arg)


def seq_file_to_counts(infile_name, fasta_file_name, count_file_name, library_fasta,
                       start_index=None, stop_index=None, test=False):
    print_now(f'Processing {infile_name}')

    file_type = None

    for file_match_string, file_type_abbr in ACCEPTED_FILE_TYPES.items():
        if fnmatch.fnmatch(infile_name, file_match_string):
            file_type = file_type_abbr
            break

    if file_type == 'fqgz':
        infile = gzip.open(infile_name)
    elif file_type == 'fq':
        infile = open(infile_name)
    elif file_type == 'fa':
        infile = open(infile_name)
    else:
        raise ValueError('Sequencing file type not recognized!')

    seq_to_id_dict, ids_to_readcount_dict, expected_read_length = parse_library_fasta(library_fasta)

    cur_read = 0
    num_aligning = 0

    with open(fasta_file_name, 'w') as unalignedFile:
        for i, fastqLine in enumerate(infile):
            if i % 4 != 1:
                continue

            else:
                seq = fastqLine.strip()[start_index:stop_index]

                if i == 1 and len(seq) != expected_read_length:
                    raise ValueError('Trimmed read length does not match expected reference read length')

                if seq in seq_to_id_dict:
                    for seqId in seq_to_id_dict[seq]:
                        ids_to_readcount_dict[seqId] += 1

                    num_aligning += 1

                else:
                    unalignedFile.write(f'> {i}\n{seq}\n')

                cur_read += 1

                # allow test runs using only the first N reads from the fastq file
                if test and cur_read >= TEST_LINES:
                    break

    with open(count_file_name, 'w') as countFile:
        for countTup in (sorted(zip(list(ids_to_readcount_dict.keys()), list(ids_to_readcount_dict.values())))):
            countFile.write(f'{countTup[0]}\t{countTup[1]}\n')

    print_now(f'Done processing {infile_name}')
    print_now(f'Counts file: {count_file_name}')

    return cur_read, num_aligning, num_aligning * 100.0 / cur_read


def parse_library_fasta(library_fasta):
    """Map File to Counts File Functions """
    seq_to_ids = dict()
    ids_to_readcounts = dict()
    read_lengths = []

    cur_seq_id = None
    cur_seq = None

    with open(library_fasta) as infile:
        for line in infile:
            if line[0] == '>':
                cur_seq_id = line.strip()[1:]
                cur_seq = ''
            else:
                cur_seq += line.strip().upper()

                if cur_seq not in seq_to_ids:
                    seq_to_ids[cur_seq] = []
                seq_to_ids[cur_seq].append(cur_seq_id)
                ids_to_readcounts[cur_seq_id] = 0
                read_lengths.append(len(cur_seq))

    if len(seq_to_ids) == 0 or len(ids_to_readcounts) == 0 or read_lengths[0] == 0:
        raise ValueError('library fasta could not be parsed or contains no sequences')
    elif max(read_lengths) != min(read_lengths):
        print(min(read_lengths), max(read_lengths))
        raise ValueError('library reference sequences are of inconsistent lengths')

    return seq_to_ids, ids_to_readcounts, read_lengths[0]


# Utility Functions

def parse_seq_file_names(file_name_list):
    infile_list = []
    outfile_base_list = []
    if isinstance(file_name_list, str):
        file_name_list = [file_name_list]

    for inputFileName in file_name_list:  # iterate through entered filenames for sequence files
        for file_name in glob.glob(inputFileName):  # generate all possible files given wildcards
            for match_string in ACCEPTED_FILE_TYPES:  # iterate through allowed filetypes
                if fnmatch.fnmatch(file_name, match_string):
                    infile_list.append(file_name)
                    outfile_base_list.append(os.path.split(file_name)[-1].split('.')[0])

    return infile_list, outfile_base_list


def make_directory(path):
    try:
        os.makedirs(path)
    except OSError:
        # printNow(path + ' already exists')
        pass


def print_now(print_input):
    print(print_input)
    sys.stdout.flush()


if __name__ == '__main__':
    args = parse_arguments()

    infileList, outfileBaseList = parse_seq_file_names(args.Seq_File_Names)
    if len(infileList) == 0:
        sys.exit('Input error: no sequencing files found')

    try:
        seqToIdDict, idsToReadcountDict, expectedReadLength = parse_library_fasta(args.Library_Fasta)

        print_now(f'Library file loaded successfully:\n\t'
                  '{:,} elements ({:,} unique sequences)\t{}bp reads expected'
                  .format(len(idsToReadcountDict), len(seqToIdDict), expectedReadLength))

    except IOError:
        sys.exit('Input error: library fasta file not found')

    except ValueError as err:
        sys.exit('Input error: ' + err.args[0])

    trimmedFastaPath = os.path.join(args.Out_File_Path, 'unaligned_reads')
    make_directory(trimmedFastaPath)
    countFilePath = os.path.join(args.Out_File_Path, 'count_files')
    make_directory(countFilePath)

    fastaFileNameList = [outfileName + '_unaligned.fa' for outfileName in outfileBaseList]
    fastaFilePathList = [os.path.join(trimmedFastaPath, fastaFileName) for fastaFileName in fastaFileNameList]
    countFilePathList = [
        os.path.join(countFilePath, outfileName + '_' + os.path.split(args.Library_Fasta)[-1] + '.counts') for
        outfileName in outfileBaseList]

    numProcessors = max(args.processors, 1)
    pool = multiprocessing.Pool(min(len(infileList), numProcessors))

    try:
        resultList = parallel_seq_file_to_counts_parallel(infileList, fastaFilePathList,
                                                          countFilePathList, pool, args.Library_Fasta, args.trim_start,
                                                          args.trim_end, args.test)
    except ValueError as err:
        sys.exit('Error while processing sequencing files: ' + ' '.join(err.args))

    for filename, result in resultList:
        print(filename + f"\n\t{result[0]:,} reads\t{result[1]:,} aligning ({result[2]:.2f}%)")

    pool.close()
    pool.join()

    print_now('Done processing all sequencing files')
