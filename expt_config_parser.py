from configparser import ConfigParser
import os


# Parse and validate the input of an experiment config file
# output a dict with all of the parameters needed to process experiments
def parse_expt_config(config_file, libraries_to_sublibraries_dict):
    parser = ConfigParser()
    results = parser.read(config_file)
    if len(results) == 0:
        return None, 1, 'Experiment config file not found'

    # output variables
    param_dict = dict()
    exit_status = 0
    warning_string = ''

    # check all sections
    expected_sections = {'experiment_settings', 'library_settings', 'counts_files',
                         'filter_settings', 'sgrna_analysis',
                         'growth_values', 'gene_analysis'}

    parsed_sections = set(parser.sections())

    if len(expected_sections) != len(parsed_sections) and len(expected_sections) != len(
            expected_sections.intersection(parsed_sections)):
        return param_dict, 1, \
               'Config file does not have all required sections or has extraneous sections!\nExpected:' + ','.join(
                expected_sections) + '\nFound:' + ','.join(parsed_sections)

    # experiment settings
    if parser.has_option('experiment_settings', 'output_folder'):
        param_dict['output_folder'] = parser.get('experiment_settings',
                                                 'output_folder')
        # TODO: ways to check this is a valid path?
    else:
        warning_string += 'No output folder specified, defaulting to current directory\n.'
        param_dict['output_folder'] = os.getcwd()

    if parser.has_option('experiment_settings', 'experiment_name'):
        param_dict['experiment_name'] = parser.get('experiment_settings', 'experiment_name')
    else:
        warning_string += 'No experiment name specified, defaulting to \'placeholder_expt_name\'\n.'
        param_dict['experiment_name'] = 'placeholder_expt_name'

    # library settings
    library_dict = libraries_to_sublibraries_dict
    if parser.has_option('library_settings', 'library'):
        parsed_library = parser.get('library_settings', 'library')

        if parsed_library.lower() in library_dict:
            param_dict['library'] = parsed_library.lower()
        else:
            warning_string += f'Library name \"{parsed_library}\" not recognized\n'
            exit_status += 1

    else:
        warning_string += 'No library specified\n'
        exit_status += 1

    if 'library' in param_dict:
        if parser.has_option('library_settings', 'sublibraries'):
            parsed_sub_list = parser.get('library_settings', 'sublibraries').strip().split('\n')

            param_dict['sublibraries'] = []

            for sub in parsed_sub_list:
                sub = sub.lower()
                if sub in library_dict[param_dict['library']]:
                    param_dict['sublibraries'].append(sub)

                else:
                    warning_string += f'Sublibrary {sub} not recognized\n'

        else:
            param_dict['sublibraries'] = library_dict[param_dict['library']]

    # counts files
    if parser.has_option('counts_files', 'counts_file_string'):
        counts_file_string = parser.get('counts_files', 'counts_file_string').strip()

        param_dict['counts_file_list'] = []

        for string_line in counts_file_string.split('\n'):
            string_line = string_line.strip()

            if len(string_line.split(':')) != 2 or len(string_line.split('|')) != 2:
                warning_string += 'counts file entry could not be parsed: ' + string_line + '\n'
                exit_status += 1

            else:
                parsed_path = string_line.split(':')[0]

                if not os.path.isfile(parsed_path):
                    warning_string += 'Counts file not found: ' + parsed_path + '\n'
                    exit_status += 1

                condition, replicate = string_line.split(':')[1].split('|')

                param_dict['counts_file_list'].append((condition, replicate, parsed_path))

    else:
        warning_string += 'No counts files entered\n'
        exit_status += 1

    # filter settings
    filter_options = ['either', 'both']
    if parser.has_option('filter_settings', 'filter_type') and parser.get('filter_settings',
                                                                          'filter_type').lower() in filter_options:
        param_dict['filter_type'] = parser.get('filter_settings', 'filter_type').lower()
    else:
        warning_string += 'Filter type not set or not recognized, defaulting to \'either\'\n'
        param_dict['filter_type'] = 'either'

    if parser.has_option('filter_settings', 'minimum_reads'):
        try:
            param_dict['minimum_reads'] = parser.getint('filter_settings', 'minimum_reads')
        except ValueError:
            # recommended value is 50 but seems arbitrary to default to that
            warning_string += 'Minimum read value not an integer, defaulting to 0\n'
            param_dict['minimum_reads'] = 0
    else:
        # recommended value is 50 but seems arbitrary to default to that
        warning_string += 'Minimum read value not found, defaulting to 0\n'
        param_dict['minimum_reads'] = 0

    # sgRNA Analysis
    if parser.has_option('sgrna_analysis', 'condition_string'):
        condition_string = parser.get('sgrna_analysis', 'condition_string').strip()

        param_dict['condition_tuples'] = []

        if 'counts_file_list' in param_dict:
            expected_conditions = set(zip(*param_dict['counts_file_list'])[0])
        else:
            expected_conditions = []

        entered_conditions = set()

        for condition_string_line in condition_string.split('\n'):
            condition_string_line = condition_string_line.strip()

            if len(condition_string_line.split(':')) != 3:
                warning_string += 'Phenotype condition line not understood: ' + condition_string_line + '\n'
                exit_status += 1
            else:
                phenotype, condition1, condition2 = condition_string_line.split(':')

                if condition1 not in expected_conditions or condition2 not in expected_conditions:
                    warning_string += 'One of the conditions entered does not correspond to a counts file: ' + \
                                      condition_string_line + '\n'
                    exit_status += 1
                else:
                    param_dict['condition_tuples'].append((phenotype, condition1, condition2))
                    entered_conditions.add(condition1)
                    entered_conditions.add(condition2)

        if len(param_dict['condition_tuples']) == 0:
            warning_string += 'No phenotype score/condition pairs found\n'
            exit_status += 1

        unused_conditions = list(expected_conditions - entered_conditions)
        if len(unused_conditions) > 0:
            warning_string += 'Some conditions assigned to counts files will not be incorporated in sgRNA analysis:\n' \
                             + ','.join(unused_conditions) + '\n'

    else:
        warning_string += 'No phenotype score/condition pairs entered\n'
        exit_status += 1

    pseudocount_options = ['zeros only', 'all values', 'filter out']
    if parser.has_option('sgrna_analysis', 'pseudocount_behavior') and \
            parser.get('sgrna_analysis', 'pseudocount_behavior').lower() in pseudocount_options:
        param_dict['pseudocount_behavior'] = parser.get('sgrna_analysis', 'pseudocount_behavior').lower()
    else:
        warning_string += 'Pseudocount behavior not set or not recognized, defaulting to \'zeros only\'\n'
        param_dict['pseudocount_behavior'] = 'zeros only'

    if parser.has_option('sgrna_analysis', 'pseudocount'):
        try:
            param_dict['pseudocount'] = parser.getfloat('sgrna_analysis', 'pseudocount')
        except ValueError:
            warning_string += 'Pseudocount value not an number, defaulting to 0.1\n'
            param_dict['pseudocount'] = 0.1
    else:
        warning_string += 'Pseudocount value not found, defaulting to 0.1\n'
        param_dict['pseudocount'] = 0.1

    # Growth Values
    if parser.has_option('growth_values', 'growth_value_string') and len(
            parser.get('growth_values', 'growth_value_string').strip()) != 0:
        growth_value_string = parser.get('growth_values', 'growth_value_string').strip()

        if 'condition_tuples' in param_dict and 'counts_file_list' in param_dict:
            expected_comparisons = set(zip(*param_dict['condition_tuples'])[0])
            expected_replicates = set(zip(*param_dict['counts_file_list'])[1])

            expected_tuple_list = []

            for comp in expected_comparisons:
                for rep in expected_replicates:
                    expected_tuple_list.append((comp, rep))
        else:
            expected_tuple_list = []

        entered_tuple_list = []
        growth_value_tuples = []

        for growth_value_line in growth_value_string.split('\n'):
            growth_value_line = growth_value_line.strip()

            linesplit = growth_value_line.split(':')

            if len(linesplit) != 3:
                warning_string += 'Growth value line not understood: ' + growth_value_line + '\n'
                exit_status += 1
                continue

            comparison = linesplit[0]
            replicate = linesplit[1]

            try:
                growth_val = float(linesplit[2])
            except ValueError:
                warning_string += 'Growth value not a number: ' + growth_value_line + '\n'
                exit_status += 1
                continue

            cur_tup = (comparison, replicate)
            if cur_tup in expected_tuple_list:
                if cur_tup not in entered_tuple_list:
                    entered_tuple_list.append(cur_tup)
                    growth_value_tuples.append((comparison, replicate, growth_val))

                else:
                    warning_string += ':'.join(cur_tup) + ' has multiple growth values entered\n'
                    exit_status += 1
            else:
                warning_string += ':'.join(
                    cur_tup) + ' was not expected given the specified counts file assignments and sgRNA phenotypes\n'
                exit_status += 1

        # because we enforced no duplicates or unexpected values these should match up unless there
        # were values not entered require all growth values to be explictly entered if some were
        if len(entered_tuple_list) != len(expected_tuple_list):
            warning_string += 'Growth values were not entered for all expected comparisons/replicates. Expected: ' + \
                             ','.join([':'.join(tup) for tup in expected_tuple_list]) + '\nEntered: ' + \
                             ','.join([':'.join(tup) for tup in entered_tuple_list]) + '\n'
            exit_status += 1
        else:
            param_dict['growth_value_tuples'] = growth_value_tuples

    else:
        warning_string += 'No growth values--all phenotypes will be reported as log2enrichments\n'

        param_dict['growth_value_tuples'] = []

        if 'condition_tuples' in param_dict and 'counts_file_list' in param_dict:
            expected_comparisons = set(zip(*param_dict['condition_tuples'])[0])
            expected_replicates = set(zip(*param_dict['counts_file_list'])[1])

            for comp in expected_comparisons:
                for rep in expected_replicates:
                    param_dict['growth_value_tuples'].append((comp, rep, 1))

    # Gene Analysis
    if parser.has_option('gene_analysis', 'collapse_to_transcripts'):
        try:
            param_dict['collapse_to_transcripts'] = parser.getboolean('gene_analysis', 'collapse_to_transcripts')
        except ValueError:
            warning_string += 'Collapse to transcripts entry not a recognized boolean value\n'
            exit_status += 1
    else:
        param_dict['collapse_to_transcripts'] = True
        warning_string += 'Collapse to transcripts defaulting to True\n'

    # pseudogene parameters
    if parser.has_option('gene_analysis', 'generate_pseudogene_dist'):
        param_dict['generate_pseudogene_dist'] = parser.get('gene_analysis', 'generate_pseudogene_dist').lower()

        if param_dict['generate_pseudogene_dist'] not in ['auto', 'manual', 'off']:
            warning_string += 'Generate pseudogene dist entry not a recognized option\n'
            exit_status += 1
    else:
        param_dict['generate_pseudogene_dist'] = False
        warning_string += 'Generate pseudogene dist defaulting to False\n'

    if 'generate_pseudogene_dist' in param_dict and param_dict['generate_pseudogene_dist'] == 'manual':
        if parser.has_option('gene_analysis', 'pseudogene_size'):
            try:
                param_dict['pseudogene_size'] = parser.getint('gene_analysis', 'pseudogene_size')
            except ValueError:
                warning_string += 'Pseudogene size entry not a recognized integer value\n'
                exit_status += 1
        else:
            warning_string += 'No pseudogene size provided\n'
            exit_status += 1

        if parser.has_option('gene_analysis', 'num_pseudogenes'):
            try:
                param_dict['num_pseudogenes'] = parser.getint('gene_analysis', 'num_pseudogenes')
            except ValueError:
                warning_string += 'Pseudogene number entry not a recognized integer value\n'
                exit_status += 1
        else:
            warning_string += 'No pseudogene size provided\n'

    # list possible analyses in param dict as dictionary with keys = analysis and values = analysis-specific params

    param_dict['analyses'] = dict()

    # analyze by average of best n
    if parser.has_option('gene_analysis', 'calculate_ave'):
        try:
            if parser.getboolean('gene_analysis', 'calculate_ave'):
                param_dict['analyses']['calculate_ave'] = []
        except ValueError:
            warning_string += 'Calculate ave entry not a recognized boolean value\n'
            exit_status += 1

        if 'calculate_ave' in param_dict['analyses']:
            if parser.has_option('gene_analysis', 'best_n'):
                try:
                    param_dict['analyses']['calculate_ave'].append(parser.getint('gene_analysis', 'best_n'))
                except ValueError:
                    warning_string += 'Best_n entry not a recognized integer value\n'
                    exit_status += 1
            else:
                warning_string += 'No best_n value provided for average analysis function\n'
                exit_status += 1
    else:
        warning_string += 'Best n average analysis not specified, defaulting to False\n'

    # analyze by Mann-Whitney
    if parser.has_option('gene_analysis', 'calculate_mw'):
        try:
            if parser.getboolean('gene_analysis', 'calculate_mw'):
                param_dict['analyses']['calculate_mw'] = []
        except ValueError:
            warning_string += 'Calculate Mann-Whitney entry not a recognized boolean value\n'
            exit_status += 1

    # analyze by K-S, skipping for now

    # analyze by nth best sgRNA
    if parser.has_option('gene_analysis', 'calculate_nth'):
        try:
            if parser.getboolean('gene_analysis', 'calculate_nth'):
                param_dict['analyses']['calculate_nth'] = []
        except ValueError:
            warning_string += 'Calculate best Nth sgRNA entry not a recognized boolean value\n'
            exit_status += 1

        if 'calculate_nth' in param_dict['analyses']:
            if parser.has_option('gene_analysis', 'nth'):
                try:
                    param_dict['analyses']['calculate_nth'].append(parser.getint('gene_analysis', 'nth'))
                except ValueError:
                    warning_string += 'Nth best sgRNA entry not a recognized integer value\n'
                    exit_status += 1
            else:
                warning_string += 'No Nth best value provided for that analysis function\n'
                exit_status += 1
    else:
        warning_string += 'Nth best sgRNA analysis not specified, defaulting to False\n'

    if len(param_dict['analyses']) == 0:
        warning_string += 'No analyses selected to compute gene scores\n'  # should this raise exit_status?

    return param_dict, exit_status, warning_string


def parse_library_config(lib_config_file):
    """Parse the library configuration file to get the available libraries, sublibraries,
     and corresponding library table files"""

    parser = ConfigParser()
    result = parser.read(lib_config_file)
    if len(result) == 0:
        raise ValueError('Library config file not found')

    libraries_to_sublibraries = dict()
    libraries_to_tables = dict()
    for library in parser.sections():
        table_file = parser.get(library, 'filename').strip()
        libraries_to_tables[library.lower()] = table_file

        sublibrary_list = parser.get(library, 'sublibraries').strip().split('\n')
        libraries_to_sublibraries[library.lower()] = [sub.strip().lower() for sub in sublibrary_list]

    if len(libraries_to_tables) == 0:
        raise ValueError('Library config file empty')

    return libraries_to_sublibraries, libraries_to_tables
