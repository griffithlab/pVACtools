import unittest
import os
import py_compile

from pvactools.lib.run_argument_utils import *
from pvactools.lib.run_utils import valid_tiers
from tests.utils import *

#python -m unittest tests/test_run_utils.py
class RunUtilsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin
        cls.utils_path = os.path.join(pvactools_directory(), "pvactools", "lib", "run_argument_utils.py")

    def test_module_compiles(self):
        self.assertTrue(py_compile.compile(self.utils_path))

    def test_pvacsplice_anchors_checker(self):
        checker = pvacsplice_anchors()

        self.assertEqual(
            checker("A,D,NDA"),
            ["A", "D", "NDA"]
        )

        with self.assertRaises(Exception) as context:
            checker("Test,A")
        self.assertEqual("List element must be one of 'A', 'D', 'NDA', not Test", str(context.exception))

    def test_float_range_checker(self):
        checker = float_range(0.0, 100.0)

        self.assertEqual(
            checker("0.5"),
            0.5
        )

        with self.assertRaises(Exception) as context:
            checker("Test")
        self.assertEqual("must be a floating point number", str(context.exception))

        with self.assertRaises(Exception) as context:
            checker("102.0")
        self.assertEqual("must be in range [0.0 .. 100.0]", str(context.exception))

    def test_aggregate_report_evaluation_checker(self):
        checker = aggregate_report_evaluations()

        self.assertEqual(
            checker("Accept,Pending"),
            ["Accept", "Pending"]
        )

        with self.assertRaises(Exception) as context:
            checker("Test,Accept")
        self.assertEqual("Invalid evaluation 'Test'. Valid values are: Accept, Reject, Pending, Review", str(context.exception))

    def test_transcript_prioritization_strategy_checker(self):
        checker = transcript_prioritization_strategy()

        self.assertEqual(
            checker("canonical,mane_select"),
            ["canonical", "mane_select"]
        )

        with self.assertRaises(Exception) as context:
            checker("Test,canonical")
        self.assertEqual("List element must be one of 'canonical', 'mane_select', 'tsl', not Test", str(context.exception))

    def test_top_score_metric2_checker(self):
        checker = top_score_metric2()

        self.assertEqual(
            checker("ic50,combined_percentile"),
            ["ic50", "combined_percentile"]
        )

        with self.assertRaises(Exception) as context:
            checker("Test,canonical")
        self.assertEqual("List element must be one of 'ic50', 'combined_percentile', 'binding_percentile', 'immunogenicity_percentile', 'presentation_percentile', not Test", str(context.exception))

    def test_tiers_checker(self):
        checker = tiers('pvacseq')

        #GENERIC
        #test that splitting into list works
        self.assertEqual(
            checker("Pass,PoorBinder"),
            ["Pass", "PoorBinder"]
        )

        #test that nonsense word fails
        with self.assertRaises(Exception) as context:
            checker("Test,Pass")
        tiers_string = ", ".join(['"{}"'.format(x) for x in valid_tiers('pvacseq')])
        self.assertEqual(f'List element must be one of {tiers_string}, not Test', str(context.exception))

        #pvacseq
        #test that all valid values work
        checker = tiers('pvacseq')
        for tier in valid_tiers('pvacseq'):
            self.assertEqual(checker(tier), [tier])
        #check that tiers unique to other pipelines fail
        with self.assertRaises(Exception) as context:
            checker("LowReadSupport")

        #pvacsplice
        #test that all valid values work
        checker = tiers('pvacsplice')
        for tier in valid_tiers('pvacsplice'):
            self.assertEqual(checker(tier), [tier])
        #check that tiers unique to other pipelines fail
        with self.assertRaises(Exception) as context:
            checker("LowReadSupport")

        #pvacfuse
        #test that all valid values work
        checker = tiers('pvacfuse')
        for tier in valid_tiers('pvacfuse'):
            self.assertEqual(checker(tier), [tier])
        #check that tiers unique to other pipelines fail
        with self.assertRaises(Exception) as context:
            checker("Subclonal")

        #pvacbind
        #test that all valid values work
        checker = tiers('pvacbind')
        for tier in valid_tiers('pvacbind'):
            self.assertEqual(checker(tier), [tier])
        #check that tiers unique to other pipelines fail
        with self.assertRaises(Exception) as context:
            checker("LowReadSupport")

    def downstream_sequence_length_checker(self):
        checker = downstream_sequence_length

        self.assertEqual(
            checker("10"),
            10
        )

        self.assertEqual(
            checker("full"),
            None
        )

        with self.assertRaises(Exception) as context:
            checker("Test")
        self.assertEqual("Argument needs to be a positive integer or 'full'", str(context.exception))

        with self.assertRaises(Exception) as context:
            checker("-1")
        self.assertEqual("Argument needs to be a positive integer or 'full'", str(context.exception))
