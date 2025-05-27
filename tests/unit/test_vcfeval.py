"""
Unit tests for the Haplo.vcfeval module.
"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, PropertyMock, patch

sys.path.insert(
    0,
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "src",
    ),
)

from hap_py.haplo import vcfeval


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

        # Mock arguments with actual RTG path
        self.args = MagicMock(spec=object)
        # Use the actual RTG path that findVCFEval returns
        self.args.engine_vcfeval = vcfeval.findVCFEval()
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
        with patch("hap_py.haplo.vcfeval.has_vcfeval", False), patch(
            "os.path.isfile", return_value=False
        ), patch("os.access", return_value=False):
            result = vcfeval.findVCFEval()
            self.assertEqual(result, "rtg")

        # Test when external RTG tools are found (current setup)
        # Don't mock anything to use actual path detection
        result = vcfeval.findVCFEval()
        # Should return actual path to RTG tools or "rtg" as fallback
        self.assertTrue(isinstance(result, str))
        # If RTG is found, it should be an absolute path
        if result != "rtg":
            self.assertTrue(os.path.isabs(result))
            self.assertTrue(result.endswith("rtg"))

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

        # Create a temporary directory for test output
        temp_dir = Path(tempfile.mkdtemp())
        temp_file_path = str(temp_dir / "vcfeval.result")

        # Create a proper mock for the context manager
        mock_context_manager = MagicMock()
        mock_context_manager.__enter__ = MagicMock()
        mock_context_manager.__exit__ = MagicMock(return_value=None)
        mock_context_manager.__enter__.return_value.name = temp_file_path

        # Mock the tempfile.NamedTemporaryFile
        with patch("tempfile.NamedTemporaryFile", return_value=mock_context_manager):
            # Should use defaults for missing parameters
            result = vcfeval.runVCFEval(
                self.test_vcf1, self.test_vcf2, self.test_output, args
            )

        shutil.rmtree(temp_dir)

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
        # Mock RTG executable existence and all other file existence checks
        with patch("os.path.isfile") as mock_isfile, patch("os.access") as mock_access:
            mock_isfile.return_value = True
            mock_access.return_value = True
            mock_exists.return_value = True
            mock_isdir.return_value = True

            # For this test skip SDF creation as we're directly testing the subprocess error handling
            self.args.engine_vcfeval_template = self.test_template

            # Create a proper temp file path
            temp_file_path = os.path.join(self.temp_dir, "vcfeval.result_mock")

            # Create a proper mock for the context manager
            mock_context_manager = MagicMock()
            mock_context_manager.__enter__ = MagicMock()
            mock_context_manager.__exit__ = MagicMock(return_value=None)
            mock_context_manager.__enter__.return_value.name = temp_file_path

            with patch(
                "tempfile.NamedTemporaryFile", return_value=mock_context_manager
            ):
                # Test vcfeval command failure - mock a failing subprocess
                process_mock = MagicMock()
                process_mock.returncode = 1
                process_mock.communicate.return_value = ("", "Error in command")
                mock_popen.return_value = process_mock

                # Create mock output directory so file checks don't interfere
                mock_out_path = os.path.join(self.temp_dir, "vcfeval.result_mock")
                os.makedirs(mock_out_path, exist_ok=True)

                # The subprocess should fail and raise SubprocessError
                with self.assertRaises(subprocess.SubprocessError):
                    vcfeval.runVCFEval(
                        self.test_vcf1, self.test_vcf2, self.test_output, self.args
                    )

    @patch("os.path.exists")
    @patch("shutil.copy")
    @patch("subprocess.Popen")
    def test_runVCFEval_missing_output(self, mock_popen, mock_copy, mock_exists):
        """Test handling of missing output files in runVCFEval."""
        # For this test, we need to mock that the RTG executable exists and is executable
        with patch("os.path.isfile") as mock_isfile, patch("os.access") as mock_access:
            # Mock RTG executable existence checks
            def mock_isfile_side_effect(path):
                if "rtg" in path and path == self.args.engine_vcfeval:
                    return True
                return True  # Mock other file existence checks as True

            mock_isfile.side_effect = mock_isfile_side_effect
            mock_access.return_value = True

            # Setup mocks - all files exist except the final output
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
