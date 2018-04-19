import unittest

from fastqgz_to_counts import parse_library_fasta, reverse_compliment


class TestFastqgzToCounts(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.fasta_file_test = '../library_reference/CRISPRi_v2_human.trim_1_29_forward.fa'

    def test_find_start_no_max_exon_len(self):
        # print(os.getcwd())
        seq_to_id_dict, ids_to_readcount_dict, expected_read_length = parse_library_fasta(self.fasta_file_test)

        # first element:
        self.assertIn('GAGACCCAGCGCTAACCAGGTTTAAGAG', seq_to_id_dict)
        self.assertEqual(seq_to_id_dict['GAGACCCAGCGCTAACCAGGTTTAAGAG'], ['A1BG_-_58858617.23-P1'])

        self.assertEqual(len(seq_to_id_dict), 205648)
        self.assertEqual(expected_read_length, 28)

        self.assertEqual(len(ids_to_readcount_dict), 209070)

        # last element:
        self.assertIn('AGAAGCCAACCTCGCGTTAGTTTAAGAG', seq_to_id_dict)
        self.assertEqual(seq_to_id_dict['AGAAGCCAACCTCGCGTTAGTTTAAGAG'], ['non-targeting_03789'])

    def test_reverse_compliment(self):

        self.assertEqual('ACGT', reverse_compliment('ACGT'))

        self.assertEqual('GCTACGT', reverse_compliment('ACGTAGC'))

        self.assertEqual('GAGACGT', reverse_compliment('acgtctc'))

        self.assertEqual('AAGCCCGTAGCA', reverse_compliment('tGCtACGggctT'))

        self.assertEqual('ATGCTGCA', reverse_compliment('tgcaGCAT'))

        st = b'TGCCATTTATGTGAGATAAGTTTAAGAG'.decode('utf-8')
        self.assertEqual('CTCTTAAACTTATCTCACATAAATGGCA', reverse_compliment(st))
