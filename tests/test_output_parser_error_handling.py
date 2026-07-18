import unittest

from pvactools.lib.output_parser import DefaultOutputParser


class OutputParserErrorHandlingTests(unittest.TestCase):
    def test_missing_missense_wildtype_match_raises_contextual_error(self):
        parser = DefaultOutputParser(
            input_iedb_files=[],
            input_tsv_file=None,
            key_file=None,
            output_file=None,
            sample_name='test',
        )
        result = {
            'allele': 'HLA-A*02:01',
            'mt_epitope_seq': 'ABCDEFGHI',
            'tsv_index': '1.GENE.TRANSCRIPT.missense.10A/T',
        }
        wt_results = {
            '2': {'wt_epitope_seq': 'ABXDEFGHI'},
            '3': {'wt_epitope_seq': 'ABCXEFGHI'},
        }

        with self.assertRaisesRegex(
            ValueError,
            (
                r'position 1.*TSV index '
                r'1\.GENE\.TRANSCRIPT\.missense\.10A/T.*'
                r'HLA-A\*02:01.*ABCDEFGHI.*'
                r'Available wildtype positions: 2, 3'
            ),
        ):
            parser.match_wildtype_and_mutant_entry_for_missense(
                result,
                '1',
                wt_results,
                None,
            )


if __name__ == '__main__':
    unittest.main()
