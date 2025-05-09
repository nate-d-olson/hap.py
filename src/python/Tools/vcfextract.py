# coding=utf-8
#
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt

"""
VCF Extractor

This module provides functions for extracting specific fields from VCF files.
"""

import gzip
import json
import logging
import os
import re
import subprocess
import tempfile
import time
from typing import List, Dict, Optional, Union, Tuple


# Python 3 compatibility for file handling
def open_file(filename: str, mode: str = "r") -> Union[gzip.GzipFile, open]:
    """
    Helper function to open files in the correct mode for both text and binary.

    Args:
        filename: Path to the file to open
        mode: Mode to open the file in (default: "r")

    Returns:
        File object in the appropriate mode
    """
    if "b" in mode:
        return open(filename, mode)
    return open(filename, mode, encoding="utf-8")


def field(val: str) -> Optional[Union[int, float, List]]:
    """
    Parse a field value into its appropriate type.

    Args:
        val: String value to parse

    Returns:
        Parsed value as int, float, or list of parsed values, or None if parsing fails
    """
    if "," in val:
        return [field(v) for v in val.split(",")]
    
    try:
        return int(val)
    except (ValueError, TypeError):
        pass

    try:
        return float(val)
    except (ValueError, TypeError):
        pass

    return None


def get_info(istr: str) -> Dict[str, Union[bool, int, float, List]]:
    """
    Parse VCF INFO field into a dictionary.

    Args:
        istr: INFO field string

    Returns:
        Dictionary of INFO field key-value pairs
    """
    result = {}
    for item in istr.split(";"):
        key, *value = [q.strip() for q in item.split("=", 1)]
        result[key] = True if not value else field(value[0])
    return result


def get_formats(fstr: str, fsample: str) -> Dict[str, Union[int, float, List]]:
    """
    Parse VCF FORMAT field and sample data into a dictionary.

    Args:
        fstr: FORMAT field string
        fsample: Sample data string

    Returns:
        Dictionary mapping FORMAT fields to their parsed values
    """
    format_fields = fstr.split(":")
    sample_values = fsample.split(":")
    
    return {field: field(value) for field, value in zip(format_fields, sample_values)}


def split_index(ffield: str) -> Tuple[str, Optional[int]]:
    """
    Split a field name that includes an array index.

    Args:
        ffield: Field name potentially containing an array index (e.g., "QUAL[0]")

    Returns:
        Tuple of (field_name, index) where index is None if no index was specified
    """
    match = re.match(r"(.*)\[([0-9]+)\]", ffield)
    if match:
        return match.group(1), int(match.group(2))
    return ffield, None


def vcf_extract(vcfname: str, features: List[str], filterfun: Optional[callable] = None) -> Generator[List[Optional[Union[str, int, float, List]]], None, None]:
    """
    Extract specified fields from a VCF file.

    Args:
        vcfname: Path to the VCF file
        features: List of field names to extract
        filterfun: Optional function to filter VCF lines

    Yields:
        List of extracted field values for each VCF record
    """

    try:
        file_opener = gzip.open if vcfname.endswith(".gz") else open_file
        with file_opener(vcfname, "rt", encoding="utf-8") as file_handle:
            feature_index = [split_index(f) for f in features]

            start_time = time.time()
            last_update = start_time
            record_count = 0
            
            for line in file_handle:
                if line.startswith("#"):
                    continue
                
                line = line.rstrip("\n")
                if filterfun and filterfun(line):
                    continue

                record_count += 1
                current_time = time.time()

                if current_time - last_update > 10:
                    last_update = current_time
                    elapsed = current_time - start_time
                    
                    try:
                        records_per_second = record_count / elapsed
                        microseconds_per_record = 1000000.0 / records_per_second
                    except ZeroDivisionError:
                        records_per_second = 0
                        microseconds_per_record = 0
                    
                    logging.info(
                        f"Since start: {record_count} records in {elapsed:.2f} seconds,"
                        f" {microseconds_per_record:.2f} us/record"
                    )

                fields = line.split("\t")
                current_values = []
                current_info = None
                current_formats = {}
                
                for index, feature in enumerate(features):
                    try:
                        if feature.lower().startswith("chr"):
                            current_values.append(fields[0])
                        elif feature.lower().startswith("pos"):
                            try:
                                current_values.append(int(fields[1]))
                            except ValueError:
                                logging.warning(f"Invalid position value: {fields[1]}")
                                current_values.append(None)
                        elif feature.lower().startswith("id"):
                            current_values.append(fields[2])
                        elif feature.lower().startswith("ref"):
                            current_values.append(fields[3])
                        elif feature.lower().startswith("alt"):
                            val = fields[4].split(",")
                            if feature_index[index][1] is not None:
                                if feature_index[index][1] < len(val):
                                    val = val[feature_index[index][1]]
                                else:
                                    val = None
                            current_values.append(val)
                        elif feature.lower().startswith("qual"):
                            try:
                                current_values.append(float(fields[5]))
                            except (ValueError, TypeError) as e:
                                logging.warning(f"Error parsing QUAL value: {e}")
                                current_values.append(None)
                        elif feature.lower().startswith("fil"):
                            filters = [] if fields[6] in {"PASS", "."} else fields[6].split(",")
                            
                            if feature_index[index][1] is not None:
                                try:
                                    current_values.append(filters[feature_index[index][1]])
                                except IndexError:
                                    current_values.append(None)
                            else:
                                current_values.append(filters)
                        elif feature.startswith("I."):
                            if current_info is None:
                                try:
                                    current_info = get_info(fields[7])
                                except Exception as e:
                                    logging.error(f"Error parsing INFO field: {e}")
                                    current_info = {}
                            
                            try:
                                info_field, index = feature_index[index]
                                value = current_info[info_field[2:]]
                                
                                if index is not None:
                                    try:
                                        current_values.append(value[index])
                                    except IndexError:
                                        logging.warning(f"Index {index} out of range for value: {value}")
                                        current_values.append(None)
                                else:
                                    current_values.append(value)
                            except (KeyError, IndexError) as e:
                                logging.warning(f"Error accessing INFO field {feature}: {e}")
                                current_values.append(None)
                        elif feature.startswith("S."):
                            field_name, index = feature_index[index]
                            sample_index = int(field_name.split(".", 3)[1])
                            field = field_name.split(".", 3)[2]
                            
                            try:
                                if sample_index not in current_formats:
                                    current_formats[sample_index] = get_formats(
                                        fields[8], fields[8 + sample_index]
                                    )
                                value = current_formats[sample_index][field]
                                
                                if index is not None:
                                    try:
                                        current_values.append(value[index])
                                    except IndexError:
                                        logging.warning(f"Index {index} out of range for value: {value}")
                                        current_values.append(None)
                                else:
                                    current_values.append(value)
                            except (KeyError, IndexError) as e:
                                logging.warning(
                                    f"Error accessing sample field {field} for sample {sample_index}: {e}"
                                )
                                current_values.append(None)
                        else:
                            current_values.append(feature)
                    except Exception as e:
                        logging.error(f"Error processing feature {feature}: {e}")
                        current_values.append(None)
                
                yield current_values
    except (IOError, OSError) as e:
        logging.error(f"Error reading VCF file {vcfname}: {e}")
        raise
        feature_index = [split_index(f) for f in features]

        start_time = time.time()
        last_update = start_time
        record_count = 0
        
        for line in ff:
            if line.startswith("#"):
                continue
            
            line = line.rstrip("\n")
            if filterfun and filterfun(line):
                continue

            record_count += 1
            current_time = time.time()

            if current_time - last_update > 10:
                last_update = current_time
                elapsed = current_time - start_time
                
                try:
                    records_per_second = record_count / elapsed
                    microseconds_per_record = 1000000.0 / records_per_second
                except ZeroDivisionError:
                    records_per_second = 0
                    microseconds_per_record = 0
                
                logging.info(
                    f"Since start: {record_count} records in {elapsed:.2f} seconds,"
                    f" {microseconds_per_record:.2f} us/record"
                )

            fields = line.split("\t")
            current_values = []
            current_info = None
            current_formats = {}
            for index, feature in enumerate(features):
                if feature.lower().startswith("chr"):
                    current_values.append(fields[0])
                elif feature.lower().startswith("pos"):
                    try:
                        current_values.append(int(fields[1]))
                    except ValueError:
                        logging.warning(f"Invalid position value: {fields[1]}")
                        current_values.append(None)
                elif feature.lower().startswith("id"):
                    current_values.append(fields[2])
                elif feature.lower().startswith("ref"):
                    current_values.append(fields[3])
                elif feature.lower().startswith("alt"):
                    val = fields[4].split(",")
                    if feature_index[index][1] is not None:
                        if feature_index[index][1] < len(val):
                            val = val[feature_index[index][1]]
                        else:
                            val = None
                    current_values.append(val)
                elif feature.lower().startswith("qual"):
                    try:
                        current_values.append(float(fields[5]))
                    except (ValueError, TypeError) as e:
                        logging.warning(f"Error parsing QUAL value: {e}")
                        current_values.append(None)
                elif feature.lower().startswith("fil"):
                    filters = [] if fields[6] in {"PASS", "."} else fields[6].split(",")
                    
                    if feature_index[index][1] is not None:
                        try:
                            current_values.append(filters[feature_index[index][1]])
                        except IndexError:
                            current_values.append(None)
                    else:
                        current_values.append(filters)
                elif feature.startswith("I."):
                    if current_info is None:
                        try:
                            current_info = get_info(fields[7])
                        except Exception as e:
                            logging.error(f"Error parsing INFO field: {e}")
                            current_info = {}
                    
                    try:
                        info_field, index = feature_index[index]
                        value = current_info[info_field[2:]]
                        
                        if index is not None:
                            try:
                                current_values.append(value[index])
                            except IndexError:
                                logging.warning(f"Index {index} out of range for value: {value}")
                                current_values.append(None)
                        else:
                            current_values.append(value)
                    except (KeyError, IndexError) as e:
                        logging.warning(f"Error accessing INFO field {feature}: {e}")
                        current_values.append(None)
                elif feature.startswith("S."):
                    field_name, index = feature_index[index]
                    sample_index = int(field_name.split(".", 3)[1])
                    field = field_name.split(".", 3)[2]
                    
                    try:
                        if sample_index not in current_formats:
                            current_formats[sample_index] = get_formats(
                                fields[8], fields[8 + sample_index]
                            )
                        value = current_formats[sample_index][field]
                        
                        if index is not None:
                            try:
                                current_values.append(value[index])
                            except IndexError:
                                logging.warning(f"Index {index} out of range for value: {value}")
                                current_values.append(None)
                        else:
                            current_values.append(value)
                    except (KeyError, IndexError) as e:
                        logging.warning(
                            f"Error accessing sample field {field} for sample {sample_index}: {e}"
                        )
                        current_values.append(None)
                else:
                    current_values.append(feature)
            yield current_values


def extract_headers(vcfname: str) -> Generator[str, None, None]:
    """
    Extract header lines from a VCF file.

    Args:
        vcfname: Path to the VCF file

    Yields:
        Header lines from the VCF file
    """
    try:
        file_opener = gzip.open if vcfname.endswith(".gz") else open_file
        with file_opener(vcfname, "rt", encoding="utf-8") as file_handle:
            for line in file_handle:
                if line.startswith("#"):
                    yield line.rstrip("\n")
                else:
                    break
    except (IOError, OSError) as e:
        logging.error(f"Error reading VCF file {vcfname}: {e}")
        raise
        for line in file_handle:
            if line.startswith("#"):
                yield line.rstrip("\n")
            else:
                break


def extract_headers_json(vcfname: str) -> Dict[str, Any]:
    """
    Extract VCF header and convert to JSON format.

    Args:
        vcfname: Path to the VCF file

    Returns:
        Dictionary containing the VCF header information

    Raises:
        Exception: If vcfhdr2json command fails
    """
    temp_file = tempfile.NamedTemporaryFile(delete=False)
    temp_file.close()
    
    try:
        # Create temporary file for vcfhdr2json output
        temp_file = tempfile.NamedTemporaryFile(delete=False)
        temp_file.close()
        
        try:
            # Run vcfhdr2json command
            process = subprocess.Popen(
                ["vcfhdr2json", vcfname, temp_file.name],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            stdout, stderr = process.communicate()

            if process.returncode != 0:
                raise Exception(f"vcfhdr2json call failed: {stdout} / {stderr}")

            # Read and parse JSON output
            with open(temp_file.name, "r", encoding="utf-8") as json_file:
                header_data = json.load(json_file)

            # Process header data
            result = {
                "aligners": [],
                "callers": []
            }
            
            for header_field in header_data.get("fields", []):
                try:
                    key = header_field.get("key")
                    values = header_field.get("values", {})
                    
                    if key == "source":
                        caller_info = ["unknown", "unknown", ""]
                        try:
                            caller_info[0] = str(values)
                            if caller_info[0].startswith("Platypus_Version_"):
                                caller_info[1] = caller_info[0][len("Platypus_Version_") :]
                                caller_info[0] = "Platypus"
                        except (ValueError, TypeError) as e:
                            logging.warning(f"Error parsing source value: {e}")
                        result["callers"].append(caller_info)
                    
                    elif key == "source_version":
                        try:
                            caller_info[1] = str(values)
                        except (ValueError, TypeError) as e:
                            logging.warning(f"Error parsing version value: {e}")
                    
                    elif key == "cmdline":
                        try:
                            caller_info[2] = str(values)
                        except (ValueError, TypeError) as e:
                            logging.warning(f"Error parsing command line value: {e}")
                    
                    elif key == "platypusOptions":
                        try:
                            caller_info[2] = str(values)
                        except (ValueError, TypeError) as e:
                            logging.warning(f"Error parsing Platypus options: {e}")
                    
                    elif key == "octopus":
                        result["callers"].append(["octopus", "unknown", str(values)])
                    
                    elif key.startswith("GATKCommandLine"):
                        caller = "GATK"
                        try:
                            caller += "-" + values.get("ID", "")
                            version = values.get("Version", "unknown")
                            options = values.get("CommandLineOptions", "")
                            
                            if any(g in caller.lower() for g in ["haplotypecaller", "unifiedgenotyper", "mutect"]):
                                result["callers"].append([caller, version, options])
                        except (ValueError, TypeError) as e:
                            logging.warning(f"Error parsing GATK information: {e}")
                    
                    elif key.startswith("SentieonCommandLine"):
                        caller = "Sentieon"
                        try:
                            caller += "-" + values.get("ID", "")
                            version = values.get("Version", "unknown")
                            
                            if any(s in caller.lower() for s in ["haplotyper"]):
                                result["callers"].append([caller, version])
                        except (ValueError, TypeError) as e:
                            logging.warning(f"Error parsing Sentieon information: {e}")
                
                except Exception as e:
                    logging.error(f"Error processing header field {key}: {e}")
                    continue
            
            return result
            
        except (json.JSONDecodeError, OSError, IOError) as e:
            logging.error(f"Error reading JSON file: {e}")
            raise
            
    except Exception as ex:
        logging.error(f"Error processing VCF header: {ex}")
        raise
        
        return header_data

    finally:
        try:
            os.unlink(temp_file.name)
        except OSError as e:
            logging.warning(f"Error cleaning up temporary file: {e}")
            return None
