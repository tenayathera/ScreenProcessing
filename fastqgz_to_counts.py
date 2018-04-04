# pipeline to generate read counts and phenotype scores directly from gzipped sequencing data

import os
import sys
import gzip
import multiprocessing
import fnmatch
import glob
import argparse


### Sequence File to Trimmed Fasta Functions ###

def parallel_seq_file_to_counts_parallel(fastq_gz_file_name_list, fasta_file_name_list, count_file_name_list, process_pool,
                                         libraryFasta, startIndex=None, stopIndex=None, test=False):
    if len(fastq_gz_file_name_list) != len(fasta_file_name_list):
        raise ValueError('In and out file lists must be the same length')

    arglist = list(
        zip(fastq_gz_file_name_list, fasta_file_name_list, count_file_name_list, [libraryFasta] * len(fasta_file_name_list),
            [startIndex] * len(fasta_file_name_list), [stopIndex] * len(fasta_file_name_list),
            [test] * len(fasta_file_name_list)))

    reads_per_file = process_pool.map(seq_file_to_counts_wrapper, arglist)

    return list(zip(count_file_name_list, reads_per_file))


def seq_file_to_counts_wrapper(arg):
    return seq_file_to_counts(*arg)


def seq_file_to_counts(infile_name, fasta_file_name, count_file_name, library_fasta,
                       start_index=None, stop_index=None, test=False):
    print_now('Processing %s' % infile_name)

    file_type = None

    for fileTup in acceptedFileTypes:
        if fnmatch.fnmatch(infile_name, fileTup[0]):
            file_type = fileTup[1]
            break

    if file_type == 'fqgz':
        infile = gzip.open(infile_name)
    elif file_type == 'fq':
        infile = open(infile_name)
    elif file_type == 'fa':
        infile = open(infile_name)
    else:
        raise ValueError('Sequencing file type not recognized!')

    seq_to_id_dict, ids_to_readcount_dict, expected_read_length = parseLibraryFasta(library_fasta)

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
                    unalignedFile.write('>%d\n%s\n' % (i, seq))

                cur_read += 1

                # allow test runs using only the first N reads from the fastq file
                if test and cur_read >= testLines:
                    break

    with open(count_file_name, 'w') as countFile:
        for countTup in (sorted(zip(list(ids_to_readcount_dict.keys()), list(ids_to_readcount_dict.values())))):
            countFile.write('%s\t%d\n' % countTup)

    print_now('Done processing %s' % infile_name)

    return cur_read, num_aligning, num_aligning * 100.0 / cur_read


### Map File to Counts File Functions ###

def parseLibraryFasta(library_fasta):
    seq_to_ids  = dict()
    ids_to_readcounts = dict()
    read_lengths = []

    cur_seq_id = ''
    cur_seq = ''

    with open(library_fasta) as infile:
        for line in infile:
            if line[0] == '>':
                if cur_seq_id != '' and cur_seq != '':
                    if cur_seq not in seq_to_ids:
                        seq_to_ids[cur_seq] = []
                    seq_to_ids[cur_seq].append(cur_seq_id)

                    ids_to_readcounts[cur_seq_id] = 0

                    read_lengths.append(len(cur_seq))

                cur_seq_id = line.strip()[1:]
                cur_seq = ''

            else:
                cur_seq += line.strip().upper()

    if len(seq_to_ids) == 0 or len(ids_to_readcounts) == 0 or read_lengths[0] == 0:
        raise ValueError('library fasta could not be parsed or contains no sequences')
    elif max(read_lengths) != min(read_lengths):
        print(min(read_lengths), max(read_lengths))
        raise ValueError('library reference sequences are of inconsistent lengths')

    return seq_to_ids, ids_to_readcounts, read_lengths[0]


### Utility Functions ###
def parse_seq_file_names(file_name_list):
    infile_list = []
    outfile_base_list = []

    for inputFileName in file_name_list:  # iterate through entered filenames for sequence files
        for filename in glob.glob(inputFileName):  # generate all possible files given wildcards
            for fileType in zip(*acceptedFileTypes)[0]:  # iterate through allowed filetypes
                if fnmatch.fnmatch(filename, fileType):
                    infile_list.append(filename)
                    outfile_base_list.append(os.path.split(filename)[-1].split('.')[0])

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


### Global variables ###
acceptedFileTypes = [('*.fastq.gz', 'fqgz'),
                     ('*.fastq', 'fq'),
                     ('*.fq', 'fq'),
                     ('*.fa', 'fa'),
                     ('*.fasta', 'fa'),
                     ('*.fna', 'fa')]

testLines = 10000

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process raw sequencing data from screens to counts files in parallel')
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
                             % testLines)

    args = parser.parse_args()
    # printNow(args)

    ###catch input mistakes###
    numProcessors = max(args.processors, 1)

    infileList, outfileBaseList = parse_seq_file_names(args.Seq_File_Names)
    if len(infileList) == 0:
        sys.exit('Input error: no sequencing files found')

    try:
        seqToIdDict, idsToReadcountDict, expectedReadLength = parseLibraryFasta(args.Library_Fasta)

        print_now('Library file loaded successfully:\n\t%.2E elements (%.2E unique sequences)\t%dbp reads expected' \
                  % (len(idsToReadcountDict), len(seqToIdDict), expectedReadLength))

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

    pool = multiprocessing.Pool(min(len(infileList), numProcessors))

    try:
        resultList = parallel_seq_file_to_counts_parallel(infileList, fastaFilePathList, countFilePathList, pool,
                                                          args.Library_Fasta, args.trim_start, args.trim_end, args.test)
    except ValueError as err:
        sys.exit('Error while processing sequencing files: ' + ' '.join(err.args))

    for filename, result in resultList:
        print(filename + ':\n\t%.2E reads\t%.2E aligning (%.2f%%)' % result)

    pool.close()
    pool.join()

    print_now('Done processing all sequencing files')
