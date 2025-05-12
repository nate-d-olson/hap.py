# coding=utf-8
#
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/sequencing/licenses/blob/master/Simplified-BSD-License.txt

"""
Caller Information Extractor

This module provides functionality for extracting caller and aligner information
from VCF and BAM files.
"""

import json
import logging
import os
import subprocess
import tempfile
from typing import Dict, List, Tuple, Any


class CallerInfo:
    """
    Class for collecting caller and aligner information.

    Attributes:
        callers: List of tuples containing (caller_name, version, parameters_string)
        aligners: List of tuples containing (aligner_name, version, parameters_string)
    """

    def __init__(self) -> None:
        """Initialize CallerInfo instance."""
        self.callers: List[Tuple[str, str, str]] = []
        self.aligners: List[Tuple[str, str, str]] = []

    def __repr__(self) -> str:
        """Return string representation of CallerInfo instance."""
        aligners_str = ",".join("/".join(x) for x in self.aligners)
        callers_str = ",".join("/".join(x) for x in self.callers)
        return f"aligners=[{aligners_str}] callers=[{callers_str}]"

    def as_dict(self) -> Dict[str, List[Dict[str, str]]]:
        """
        Convert caller and aligner information to dictionary format.

        Returns:
            Dictionary containing aligners and callers information.
        """
        key_value_list = ["name", "version", "parameters"]
        return {
            "aligners": [dict(list(zip(key_value_list, x))) for x in self.aligners],
            "callers": [dict(list(zip(key_value_list, x))) for x in self.callers],
        }

    def add_vcf(self, vcfname: str) -> None:
        """
        Add caller versions from a VCF file.

        Args:
            vcfname: Path to the VCF file.
        """
        temp_file_name = ""
        header_data: Dict[str, Any] = {}

        try:
            with tempfile.NamedTemporaryFile(
                delete=False, mode="w", encoding="utf-8"
            ) as temp_file_obj:
                temp_file_name = temp_file_obj.name

            if not temp_file_name:
                logging.error("Failed to create a temporary file.")
                return

            process = subprocess.Popen(
                ["vcfhdr2json", vcfname, temp_file_name],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            stdout, stderr = process.communicate()

            if process.returncode != 0:
                logging.error(
                    "vcfhdr2json call failed for %s: %s / %s",
                    vcfname,
                    stdout,
                    stderr,
                )
                return

            with open(temp_file_name, "r", encoding="utf-8") as json_file:
                header_data = json.load(json_file)

        except FileNotFoundError:
            logging.error(
                "vcfhdr2json command not found. Please ensure it is in your PATH."
            )
            return
        except json.JSONDecodeError as jde:
            logging.error(
                "Error decoding JSON from vcfhdr2json output for %s: %s",
                vcfname,
                jde,
            )
            return
        except Exception as ex:
            logging.error("Error processing VCF header for %s: %s", vcfname, ex)
            return
        finally:
            if temp_file_name and os.path.exists(temp_file_name):
                try:
                    os.unlink(temp_file_name)
                except OSError as e:
                    logging.warning(
                        "Error removing temp file %s: %s", temp_file_name, e
                    )

        # name, version, params
        current_caller_info: List[str] = ["unknown", "unknown", ""]
        gatk_callers = ["haplotypecaller", "unifiedgenotyper", "mutect"]
        sentieon_callers = ["haplotyper", "tnhaplotyper2", "tnscope"]
        source_info_found = False

        if not header_data or "fields" not in header_data:
            logging.warning(
                "VCF header data for %s is missing 'fields' key or is empty.",
                vcfname,
            )
            return

        for header_field in header_data.get("fields", []):
            try:
                key = str(header_field.get("key", ""))
                values = header_field.get("values", {})

                if not key:
                    continue

                if key == "source":
                    try:
                        source_val = str(values)
                        current_caller_info[0] = source_val
                        if source_val.startswith("Platypus_Version_"):
                            current_caller_info[1] = source_val[
                                len("Platypus_Version_") :
                            ]
                            current_caller_info[0] = "Platypus"
                        source_info_found = True
                    except (ValueError, TypeError) as e:
                        logging.warning(
                            "Error parsing 'source' value (%s): %s", values, e
                        )

                elif key == "source_version":
                    try:
                        current_caller_info[1] = str(values)
                        source_info_found = True
                    except (ValueError, TypeError) as e:
                        logging.warning(
                            "Error parsing 'source_version' value (%s): %s",
                            values,
                            e,
                        )

                elif key == "cmdline":
                    try:
                        current_caller_info[2] = str(values)
                        source_info_found = True
                    except (ValueError, TypeError) as e:
                        logging.warning(
                            "Error parsing 'cmdline' value (%s): %s", values, e
                        )

                elif key == "platypusOptions":
                    try:
                        # Overwrites cmdline if both present
                        current_caller_info[2] = str(values)
                        source_info_found = True
                    except (ValueError, TypeError) as e:
                        logging.warning(
                            "Error parsing 'platypusOptions' value (%s): %s",
                            values,
                            e,
                        )

                elif key == "octopus":
                    self.callers.append(("octopus", "unknown", str(values)))

                elif isinstance(values, dict) and key.startswith("GATKCommandLine"):
                    caller_name = "GATK"
                    version = "unknown"
                    options = ""
                    try:
                        sub_id = str(values.get("ID", ""))
                        if sub_id:
                            caller_name += "-" + sub_id
                        version = str(values.get("Version", "unknown"))
                        options = str(values.get("CommandLineOptions", ""))
                        # Check if the identified GATK tool is one we track
                        if any(
                            g_caller in caller_name.lower() for g_caller in gatk_callers
                        ):
                            self.callers.append((caller_name, version, options))
                    except (ValueError, TypeError, KeyError) as e:
                        logging.warning(
                            "Error parsing GATK fields for key %s (%s): %s",
                            key,
                            values,
                            e,
                        )

                elif isinstance(values, dict) and key.startswith("SentieonCommandLine"):
                    caller_name = "Sentieon"
                    version = "unknown"
                    options = ""
                    try:
                        sub_id = str(values.get("ID", ""))
                        if sub_id:
                            caller_name += "-" + sub_id
                        version = str(values.get("Version", "unknown"))
                        options = str(values.get("CommandLineOptions", ""))
                        if any(
                            s_caller in caller_name.lower()
                            for s_caller in sentieon_callers
                        ):
                            self.callers.append((caller_name, version, options))
                    except (ValueError, TypeError, KeyError) as e:
                        logging.warning(
                            "Error parsing Sentieon fields for key %s (%s): %s",
                            key,
                            values,
                            e,
                        )

            except Exception as e:
                logging.error(
                    "General error processing VCF header field (key: %s): %s",
                    str(header_field.get("key", "N/A")),
                    e,
                )

        if source_info_found:
            self.callers.append(
                (
                    current_caller_info[0],
                    current_caller_info[1],
                    current_caller_info[2],
                )
            )

    def add_bam(self, bamfile: str) -> None:
        """
        Extract aligner information from a BAM file.

        Args:
            bamfile: Path to the BAM file.
        """
        try:
            process = subprocess.Popen(
                ["samtools", "view", "-H", bamfile],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            stdout, stderr = process.communicate()

            if process.returncode != 0:
                logging.error(
                    "samtools view -H call failed for %s: %s / %s",
                    bamfile,
                    stdout,
                    stderr,
                )
                return

            for line in stdout.splitlines():
                if not line.startswith("@PG"):
                    continue

                header_dict: Dict[str, str] = {}
                try:
                    # Fields are tab-separated, first field is @PG tag
                    fields = line.strip().split("\t")[1:]
                    for field in fields:
                        if ":" in field:
                            k, v = field.split(":", 1)
                            header_dict[k] = v
                        else:
                            logging.debug(
                                "Malformed PG tag field (no colon): '%s' in line: %s",
                                field,
                                line,
                            )
                except ValueError as ve:
                    logging.warning(
                        "Unable to parse SAM/BAM PG line: '%s'. Error: %s",
                        line,
                        ve,
                    )
                    continue

                aligner_name = "unknown"
                aligner_version = "unknown"
                aligner_cl = ""

                # Prefer PN for program name, fallback to ID
                aligner_name = header_dict.get("PN", header_dict.get("ID", "unknown"))

                if aligner_name == "unknown":
                    logging.warning(
                        "Missing Program Name (PN) and ID in BAM @PG line for %s: %s",
                        bamfile,
                        line,
                    )
                    continue

                # If ID was used for name and contains a common suffix (e.g., bwa-mem),
                # try to clean it. This is a heuristic.
                if header_dict.get("PN") is None and "-" in aligner_name:
                    parts = aligner_name.split("-", 1)
                    if parts:  # Check if split produced any parts
                        aligner_name = parts[0]

                aligner_version = header_dict.get("VN", "unknown")
                aligner_cl = header_dict.get("CL", "")

                self.aligners.append((aligner_name, aligner_version, aligner_cl))

        except FileNotFoundError:
            logging.error(
                "samtools command not found. Please ensure it is in your PATH."
            )
        except Exception as e:
            logging.error("Error processing BAM file %s: %s", bamfile, e)


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, format="%(levelname)s: %(message)s")
    ci_test = CallerInfo()

    # --- Test VCF Processing ---
    VCFHDR2JSON_MOCK_SCRIPT = "./vcfhdr2json_mock.sh"
    DUMMY_VCF_PATH = "dummy_test.vcf"

    # Create a dummy vcfhdr2json script
    # Using a standard string, not an f-string, as no placeholders are needed here.
    mock_script_content = """#!/bin/sh
# Mock vcfhdr2json script
echo '{
  "fields": [
    {
      "key": "source",
      "values": "MyTestCallerV1.2.3"
    },
    {
      "key": "cmdline",
      "values": "MyTestCaller --options"
    },
    {
      "key": "GATKCommandLine",
      "values": {
        "ID": "HaplotypeCaller",
        "Version": "4.1",
        "CommandLineOptions": "-R ref.fa -I in.bam"
      }
    }
  ]
}' > "$2"
exit 0
"""
    with open(VCFHDR2JSON_MOCK_SCRIPT, "w", encoding="utf-8") as f_script:
        f_script.write(mock_script_content)
    os.chmod(VCFHDR2JSON_MOCK_SCRIPT, 0o755)

    # Create a dummy VCF file
    dummy_vcf_content = (
        "##fileformat=VCFv4.2\n"
        "##source=MyTestCallerV1.2.3\n"
        "##GATKCommandLine=<ID=HaplotypeCaller,Version=4.1,"
        'CommandLineOptions="-R ref.fa -I in.bam">\n'
        '##cmdline="MyTestCaller --options"\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t123\t.\tA\tT\t100\tPASS\t.\n"
    )
    with open(DUMMY_VCF_PATH, "w", encoding="utf-8") as f_vcf:
        f_vcf.write(dummy_vcf_content)

    original_path = os.environ.get("PATH", "")
    mock_script_dir = os.path.abspath(os.path.dirname(VCFHDR2JSON_MOCK_SCRIPT))

    # Use a specific name for the mock during test to avoid conflict
    renamed_mock_path = os.path.join(mock_script_dir, "vcfhdr2json")

    try:
        if os.path.exists(renamed_mock_path):  # Clean up if exists from previous run
            os.remove(renamed_mock_path)
        os.rename(VCFHDR2JSON_MOCK_SCRIPT, renamed_mock_path)
        os.environ["PATH"] = f"{mock_script_dir}{os.pathsep}{original_path}"

        logging.info("Testing add_vcf with %s...", DUMMY_VCF_PATH)
        ci_test.add_vcf(DUMMY_VCF_PATH)
        logging.info("Collected Caller Info after VCF: %s", ci_test.as_dict())
    except Exception as e_main:
        logging.error("Error during __main__ VCF test: %s", e_main)
    finally:
        os.environ["PATH"] = original_path
        if os.path.exists(renamed_mock_path):
            os.remove(renamed_mock_path)
        if os.path.exists(VCFHDR2JSON_MOCK_SCRIPT):  # In case rename failed
            os.remove(VCFHDR2JSON_MOCK_SCRIPT)
        if os.path.exists(DUMMY_VCF_PATH):
            os.remove(DUMMY_VCF_PATH)

    # --- Test BAM Processing (Conceptual) ---
    logging.info("Skipping add_bam test in this example due to samtools dependency.")

    logging.info("Test run finished. Final Caller Info: %s", ci_test.as_dict())
