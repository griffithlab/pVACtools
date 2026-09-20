from contextlib import redirect_stdout
import importlib
from io import StringIO
import sys
import unittest
from unittest.mock import patch


CLI_COMMANDS = [
    ('pvactools.tools.main', 'valid_alleles'),
    ('pvactools.tools.pvacseq.main', 'run'),
    ('pvactools.tools.pvacbind.main', 'run'),
    ('pvactools.tools.pvacfuse.main', 'run'),
    ('pvactools.tools.pvacsplice.main', 'run'),
    ('pvactools.tools.pvacvector.main', 'run'),
    ('pvactools.tools.pvacview.main', 'run'),
]


class CliErrorHandlingTests(unittest.TestCase):
    def test_command_attribute_errors_are_not_hidden(self):
        for module_name, command in CLI_COMMANDS:
            with self.subTest(module=module_name):
                cli_module = importlib.import_module(module_name)
                command_module = getattr(cli_module, command)
                argv = [module_name, command]
                error = AttributeError('business failure')
                with patch.object(sys, 'argv', argv), patch.object(
                    command_module,
                    'main',
                    side_effect=error,
                ):
                    with self.assertRaises(AttributeError) as context:
                        cli_module.main()
                self.assertIs(context.exception, error)

    def test_missing_commands_report_cli_error(self):
        for module_name, _ in CLI_COMMANDS:
            with self.subTest(module=module_name):
                cli_module = importlib.import_module(module_name)
                output = StringIO()
                with patch.object(sys, 'argv', [module_name]), redirect_stdout(output):
                    with self.assertRaises(SystemExit) as context:
                        cli_module.main()
                self.assertEqual(context.exception.code, -1)
                self.assertIn('Error: No command specified', output.getvalue())


if __name__ == '__main__':
    unittest.main()
