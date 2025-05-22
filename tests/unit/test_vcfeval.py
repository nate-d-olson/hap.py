"""
Unit tests for the Haplo.vcfeval module.
"""

import os
import shutil
import subprocess

# Add src/python to path for imports during tests
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, PropertyMock, patch

sys.path.insert(
    0,
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "src",
        "python",
    ),
)

from Haplo import vcfeval


class TestVCFEval(unittest.TestCase):
    """Test cases for the vcfeval module."""

    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.test_vcf1 = os.path.join(self.temp_dir, "test1.vcf.gz")
        self.test_vcf2 = os.path.join(self.temp_dir, "test2.vcf.gz")
        self.test_output = os.path.join(self.temp_dir, "output.vcf.gz")
        self.test_ref = os.path.join(self.temp_dir, "ref.fa")
        self.test_template = os.path.join(self.temp_dir, "template")

        # Create empty files for testing
        with open(self.test_vcf1, "w") as f:
            f.write("")
        with open(self.test_vcf2, "w") as f:
            f.write("")
        with open(self.test_ref, "w") as f:
            f.write("")

        # Create directory for template
        os.makedirs(self.test_template, exist_ok=True)

        # Mock arguments
        self.args = MagicMock(spec=object)
        self.args.engine_vcfeval = "rtg"
        self.args.engine_vcfeval_template = self.test_template
        self.args.ref = self.test_ref
        self.args.scratch_prefix = self.temp_dir
        self.args.threads = 1
        self.args.pass_only = False
        self.args.roc = None

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.temp_dir)

    def test_findVCFEval(self):
        """Test the findVCFEval function."""
        # Test when has_vcfeval is False
        with patch("Haplo.vcfeval.has_vcfeval", False):
            result = vcfeval.findVCFEval()
            self.assertEqual(result, "rtg")

        # Test when has_vcfeval is True but files don't exist
        with patch("Haplo.vcfeval.has_vcfeval", True), patch(
            "os.path.isfile", return_value=False
        ):
            result = vcfeval.findVCFEval()
            self.assertEqual(result, "rtg")

    @patch("subprocess.Popen")
    @patch("shutil.copy")
    @patch("os.path.exists")
    def test_runVCFEval_input_validation(self, mock_exists, mock_copy, mock_popen):
        """Test input validation in runVCFEval."""
        # Mock for temp file
        vtf_mock = MagicMock()
        name_property = PropertyMock(
            return_value=os.path.join(self.temp_dir, "vcfeval.result_mock")
        )
        type(vtf_mock).name = name_property

        with patch("tempfile.NamedTemporaryFile", return_value=vtf_mock):
            # Test missing input file
            def exists_side_effect1(path):
                if path == self.test_vcf1:
                    return False  # First file is missing
                return True

            mock_exists.side_effect = exists_side_effect1

            # Test missing input file
            with self.assertRaises(FileNotFoundError):
                vcfeval.runVCFEval(
                    self.test_vcf1, self.test_vcf2, self.test_output, self.args
                )

            # Test missing reference file
            def exists_side_effect2(path):
                if path == self.test_ref:
                    return False  # Reference is missing
                return True

            mock_exists.side_effect = exists_side_effect2

            # Test missing reference file
            with self.assertRaises(FileNotFoundError):
                vcfeval.runVCFEval(
                    self.test_vcf1, self.test_vcf2, self.test_output, self.args
                )

    @patch("subprocess.Popen")
    @patch("shutil.copy")
    @patch("os.path.exists")
    @patch("os.path.isdir")
    def test_runVCFEval_default_params(
        self, mock_isdir, mock_exists, mock_copy, mock_popen
    ):
        """Test default parameter handling in runVCFEval."""
        # Setup mocks
        mock_exists.return_value = True
        mock_isdir.return_value = True

        # Create mock output directory and files
        mock_out_dir = os.path.join(self.temp_dir, "vcfeval.result_mock")
        os.makedirs(mock_out_dir, exist_ok=True)
        mock_out_vcf = os.path.join(mock_out_dir, "output.vcf.gz")
        mock_out_tbi = os.path.join(mock_out_dir, "output.vcf.gz.tbi")

        # Create empty output files
        with open(mock_out_vcf, "w") as f:
            f.write("")
        with open(mock_out_tbi, "w") as f:
            f.write("")

        # Mock process
        process_mock = MagicMock()
        process_mock.returncode = 0
        process_mock.communicate.return_value = ("output", "")
        mock_popen.return_value = process_mock

        # Test with missing parameters
        args = MagicMock(spec=object)
        args.ref = self.test_ref
        args.engine_vcfeval = None  # Missing engine
        args.engine_vcfeval_template = self.test_template
        args.scratch_prefix = None  # Missing scratch
        args.threads = None  # Missing threads
        args.pass_only = False
        args.roc = None

        # Create property mocks for temp file
        vtf_mock = MagicMock()
        name_property = PropertyMock(return_value=mock_out_dir)
        type(vtf_mock).name = name_property

        # Mock the tempfile.NamedTemporaryFile context manager
        with patch("tempfile.NamedTemporaryFile", return_value=vtf_mock):
            # Should use defaults for missing parameters
            result = vcfeval.runVCFEval(
                self.test_vcf1, self.test_vcf2, self.test_output, args
            )

            # Verify result is correct
            self.assertEqual(result, [self.test_output, self.test_output + ".tbi"])

            # Verify engine was set to default
            self.assertEqual(args.engine_vcfeval, vcfeval.findVCFEval())

            # Verify threads was set to default
            self.assertEqual(args.threads, 1)

            # Verify scratch_prefix was set to default
            self.assertIsNotNone(args.scratch_prefix)

    @patch("subprocess.Popen")
    @patch("shutil.copy")
    @patch("os.path.exists")
    @patch("os.path.isdir")
    def test_runVCFEval_subprocess_errors(
        self, mock_isdir, mock_exists, mock_copy, mock_popen
    ):
        """Test subprocess error handling in runVCFEval."""
        # Setup mocks for file existence
        mock_exists.return_value = True
        mock_isdir.return_value = True

        # For this test skip SDF creation as we're directly testing the subprocess error handling
        self.args.engine_vcfeval_template = self.test_template

        # Mock for temp file
        vtf_mock = MagicMock()
        name_property = PropertyMock(
            return_value=os.path.join(self.temp_dir, "vcfeval.result_mock")
        )
        type(vtf_mock).name = name_property

        with patch("tempfile.NamedTemporaryFile", return_value=vtf_mock):
            # Test format command failure - we'll skip the format step and just test the vcfeval error
            process_mock = MagicMock()
            process_mock.returncode = 1
            process_mock.communicate.return_value = ("", "Error in command")
            mock_popen.return_value = process_mock

            # Create mock output files that will be checked for
            mock_out_path = os.path.join(self.temp_dir, "vcfeval.result_mock")
            os.makedirs(mock_out_path, exist_ok=True)

            with self.assertRaises(subprocess.SubprocessError):
                vcfeval.runVCFEval(
                    self.test_vcf1, self.test_vcf2, self.test_output, self.args
                )

            # Reset mock
            mock_popen.reset_mock()

            # Test vcfeval command failure with multiple calls
            # First process succeeds (format)
            process_mock1 = MagicMock()
            process_mock1.returncode = 0
            process_mock1.communicate.return_value = ("", "")

            # Second process fails (vcfeval)
            process_mock2 = MagicMock()
            process_mock2.returncode = 1
            process_mock2.communicate.return_value = ("", "Error in vcfeval")

            # Setup popen to return different mock objects on each call
            mock_popen.side_effect = [process_mock1, process_mock2]

            with self.assertRaises(subprocess.SubprocessError):
                vcfeval.runVCFEval(
                    self.test_vcf1, self.test_vcf2, self.test_output, self.args
                )

    @patch("subprocess.Popen")
    @patch("shutil.copy")
    @patch("os.path.exists")
    def test_runVCFEval_missing_output(self, mock_exists, mock_copy, mock_popen):
        """Test handling of missing output files in runVCFEval."""
        # Setup mocks
        mock_exists.side_effect = lambda x: x != os.path.join(
            self.temp_dir, "output.vcf.gz"
        )

        process_mock = MagicMock()
        process_mock.returncode = 0
        process_mock.communicate.return_value = ("", "")
        mock_popen.return_value = process_mock

        # Test missing output file
        result = vcfeval.runVCFEval(
            self.test_vcf1, self.test_vcf2, self.test_output, self.args
        )
        self.assertIsNone(result)


if __name__ == "__main__":
    unittest.main()
