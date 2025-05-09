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

This module provides functionality for extracting caller and aligner information from VCF and BAM files.
"""

import itertools
import json
import logging
import os
import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple


class CallerInfo:
    """
    Class for collecting caller and aligner information.

    Attributes:
        callers: List of tuples containing (caller, version, parameters)
        aligners: List of tuples containing (aligner, version, parameters)
    """
    """Class for collecting caller info and version"""

    def __init__(self) -> None:
        """Initialize CallerInfo instance."""
        self.callers: List[Tuple[str, str, str]] = []
        self.aligners: List[Tuple[str, str, str]] = []

    def __repr__(self) -> str:
        """Return string representation of CallerInfo instance."""
        aligners_str = ",".join("/".join(x) for x in self.aligners)
        callers_str = ",".join("/".join(x) for x in self.callers)
        return f"aligners=[{aligners_str}] callers=[{callers_str}]"

    def as_dict(self) -> Dict:
        """
        Convert caller and aligner information to dictionary format.

        Returns:
            Dictionary containing aligners and callers information
        """
        key_value_dict = ["name", "version", "parameters"]
        return {
            "aligners": [dict(zip(key_value_dict, x)) for x in self.aligners],
            "callers": [dict(zip(key_value_dict, x)) for x in self.callers],
        }

    def add_vcf(self, vcfname: str) -> None:
        """
        Add caller versions from a VCF file.

        Args:
            vcfname: Path to the VCF file

        Raises:
            Exception: If vcfhdr2json command fails
        """
        temp_file = tempfile.NamedTemporaryFile(delete=False)
        temp_file.close()
        
        try:
            process = subprocess.Popen(
                ["vcfhdr2json", vcfname, temp_file.name],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            stdout, stderr = process.communicate()

            if process.returncode != 0:
                raise Exception(f"vcfhdr2json call failed: {stdout} / {stderr}")

            with open(temp_file.name, "r", encoding="utf-8") as json_file:
                header_data = json.load(json_file)

            self._process_vcf_header(header_data)

        except Exception as ex:
            logging.error(f"Error processing VCF header: {ex}")
        finally:
            try:
                os.unlink(temp_file.name)
            except OSError as e:
                logging.warning(f"Error removing temp file: {e}")

        caller_info = ["unknown", "unknown", ""]
        gatk_callers = ["haplotypecaller", "unifiedgenotyper", "mutect"]
        sent_callers = ["haplotyper"]
        source_found = False

        for header_field in header_data["fields"]:
            try:
                key = header_field["key"]
                values = header_field.get("values", {})

                if key == "source":
                    try:
                        caller_info[0] = str(values)
                    except (ValueError, TypeError) as e:
                        logging.warning(f"Error parsing source value: {e}")
                        continue
                    
                    if caller_info[0].startswith("Platypus_Version_"):
                        caller_info[1] = caller_info[0][len("Platypus_Version_") :]
                        caller_info[0] = "Platypus"
                    source_found = True

                elif key == "source_version":
                    try:
                        caller_info[1] = str(values)
                    except (ValueError, TypeError) as e:
                        logging.warning(f"Error parsing version value: {e}")
                        continue
                    source_found = True

                elif key == "cmdline":
                    try:
                        caller_info[2] = str(values)
                    except (ValueError, TypeError) as e:
                        logging.warning(f"Error parsing command line value: {e}")
                        continue
                    source_found = True

                elif key == "platypusOptions":
                    try:
                        caller_info[2] = str(values)
                    except (ValueError, TypeError) as e:
                        logging.warning(f"Error parsing Platypus options: {e}")
                        continue
                    source_found = True

                elif key == "octopus":
                    self.callers.append(["octopus", "unknown", str(values)])

                elif key.startswith("GATKCommandLine"):
                    caller = "GATK"
                    try:
                        caller += "-" + values["ID"]
                    except (ValueError, TypeError) as e:
                        logging.warning(f"Error parsing GATK ID: {e}")
                    else:
                        version = "unknown"
                        try:
                            version = values["Version"]
                        except (ValueError, TypeError) as e:
                            logging.warning(f"Error parsing GATK version: {e}")
                        else:
                            options = ""
                            try:
                                options = values["CommandLineOptions"]
                            except (ValueError, TypeError) as e:
                                logging.warning(f"Error parsing GATK options: {e}")
                            else:
                                if any(g in caller.lower() for g in gatk_callers):
                                    self.callers.append([caller, version, options])

                elif key.startswith("SentieonCommandLine"):
                    caller = "Sentieon"
                    try:
                        caller += "-" + values["ID"]
                    except (ValueError, TypeError) as e:
                        logging.warning(f"Error parsing Sentieon ID: {e}")
                    else:
                        version = "unknown"
                        try:
                            version = values["Version"]
                        except (ValueError, TypeError) as e:
                            logging.warning(f"Error parsing Sentieon version: {e}")
                        else:
                            if any(s in caller.lower() for s in sent_callers):
                                self.callers.append([caller, version])

            except (OSError, IOError) as e:
                logging.error(f"Error reading header field: {e}")
                continue

        if source_found:
            self.callers.append(caller_info)

    def add_bam(self, bamfile: str) -> None:
        """
        Extract aligner information from a BAM file.

        Args:
            bamfile: Path to the BAM file

        Raises:
            Exception: If samtools command fails
        """
        process = subprocess.Popen(
            ["samtools", "view", "-H", bamfile],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            raise Exception(f"Samtools call failed: {stdout} / {stderr}")

        for line in stdout.split("\n"):
            if not line.startswith("@PG"):
                continue
            
            try:
                # Split header line into key-value pairs
                fields = line.split("\t")[1:]
                header_dict = dict(field.split(":", 1) for field in fields)
            except (OSError, IOError) as e:
                logging.error(f"Error reading BAM file: {e}")
                logging.warning(f"Unable to parse SAM/BAM header line: {line}")
                continue
            
            aligner_info = ["unknown", "unknown", ""]
            try:
                aligner_info[0] = header_dict["PN"]
            except KeyError as e:
                logging.warning(f"Missing program name in BAM header: {e}")
                continue
            
            try:
                aligner_info[1] = header_dict["VN"]
            except KeyError as e:
                logging.warning(f"Missing version in BAM header: {e}")
            
            try:
                aligner_info[2] = header_dict["CL"]
            except KeyError as e:
                logging.warning(f"Missing command line in BAM header: {e}")
            
            self.aligners.append(tuple(aligner_info))

            try:
                aligner_info[0] = header_dict["ID"]
                if "-" in aligner_info[0]:
                    aligner_info[0] = aligner_info[0].split("-")[0]
            except KeyError as e:
                logging.warning(f"Missing ID in BAM header: {e}")
                continue
                continue
                pass

            self.aligners.append(cp)
