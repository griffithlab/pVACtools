import unittest
import os
import py_compile
import tempfile

from pvactools.lib.run_utils import *
from tests.utils import *

#python -m unittest tests/test_run_utils.py
class RunUtilsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin
        cls.utils_path = os.path.join(pvactools_directory(), "pvactools", "lib", "run_utils.py")

    def test_module_compiles(self):
        self.assertTrue(py_compile.compile(self.utils_path))

    def test_is_preferred_transcript_accepts_boolean_values_and_strings(self):
        for column, strategy in [
            ('Canonical', 'canonical'),
            ('MANE Select', 'mane_select'),
        ]:
            for value, expected in [
                (True, True),
                (False, False),
                ('True', True),
                ('False', False),
                ('Not Run', True),
            ]:
                with self.subTest(column=column, value=value):
                    mutation = {
                        'Canonical': False,
                        'MANE Select': False,
                    }
                    mutation[column] = value
                    self.assertEqual(
                        expected,
                        is_preferred_transcript(mutation, [strategy], 1),
                    )

    def test_is_preferred_transcript_rejects_invalid_boolean_value(self):
        for column, strategy in [
            ('Canonical', 'canonical'),
            ('MANE Select', 'mane_select'),
        ]:
            with self.subTest(column=column):
                mutation = {
                    'Canonical': False,
                    'MANE Select': False,
                }
                mutation[column] = 'yes'
                with self.assertRaises(ValueError) as context:
                    is_preferred_transcript(mutation, [strategy], 1)
                self.assertEqual(
                    "Invalid value 'yes' for {!r}. Expected True, False, or 'Not Run'.".format(column),
                    str(context.exception),
                )

    def test_is_preferred_transcript_does_not_execute_malicious_value(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            marker_file = os.path.join(tmp_dir, "executed")
            mutation = {
                'Canonical': "__import__('pathlib').Path({!r}).touch() or True".format(marker_file),
                'MANE Select': False,
            }
            with self.assertRaisesRegex(ValueError, "Invalid value"):
                is_preferred_transcript(mutation, ['canonical'], 1)
            self.assertFalse(os.path.exists(marker_file))
