#!/usr/bin/env python3
#
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt

import gzip
import json
import logging
import re
import subprocess
from typing import Any, Dict, List, Optional, Union


def field(val: str) -> Union[int, float, str, List[Any]]:
    """extract field into result, guess type"""
    if "," in val:
        val_list = [field(v) for v in val.split(",")]
        return val_list
    else:
        # Try to convert to int first, then float, otherwise leave as string
        try:
            return int(val)
        except ValueError:
            try:
                return float(val)
            except ValueError:
                return val


def getInfo(istr: str) -> Dict[str, Any]:
    """Split VCF INFO String"""
    spi = istr.split(";")
    res = {}
    for x in spi:
        ax = [q.strip() for q in x.split("=", 1)]

        if len(ax) == 1:
            res[ax[0]] = True
        elif len(ax) == 2:
            res[ax[0]] = field(ax[1])
    return res


def extract_header(
    filename: str,
    extract_columns: bool = True,
    extract_info: bool = False,
    extract_formats: bool = False,
    extract_filters: bool = False,
) -> Dict[str, Any]:
    """Extract header information from a VCF file

    Args:
        filename: Path to VCF file
        extract_columns: When true, extract column headers
        extract_info: When true, extract INFO fields
        extract_formats: When true, extract FORMAT fields
        extract_filters: When true, extract FILTER fields

    Returns:
        Dictionary with header information
    """
    result = {}

    if extract_columns:
        result["columns"] = []

    if extract_info:
        result["info"] = {}

    if extract_formats:
        result["formats"] = {}

    if extract_filters:
        result["filters"] = {}

    # For compatibility with pre.py and hap.py
    result["fields"] = []
    result["samples"] = []  # Always initialize samples field

    file_handle = None
    try:
        if filename.endswith(".gz"):
            file_handle = gzip.open(filename, "rt", encoding="utf-8")
        else:
            file_handle = open(filename, encoding="utf-8")

        for line in file_handle:
            line = line.strip()
            if line.startswith("##INFO="):
                # process INFO lines
                if not extract_info and "fields" not in result:
                    continue

                info_match = re.match(
                    r"##INFO=<ID=([^,>]+),Number=([^,>]+),Type=([^,>]+),Description=\"([^\"]+)\"",
                    line,
                )
                if not info_match:
                    logging.error("Cannot parse INFO line: %s" % line)
                else:
                    info_field = {
                        "id": info_match.group(1),
                        "number": info_match.group(2),
                        "type": info_match.group(3),
                        "description": info_match.group(4),
                    }

                    if extract_info:
                        result["info"][info_match.group(1)] = info_field

                    # Add to fields list for compatibility
                    result["fields"].append({"key": "INFO", "values": info_field})

            elif line.startswith("##FORMAT="):
                # process FORMAT lines
                if not extract_formats and "fields" not in result:
                    continue

                format_match = re.match(
                    r"##FORMAT=<ID=([^,>]+),Number=([^,>]+),Type=([^,>]+),Description=\"([^\"]+)\"",
                    line,
                )
                if not format_match:
                    logging.error("Cannot parse FORMAT line: %s" % line)
                else:
                    format_field = {
                        "id": format_match.group(1),
                        "number": format_match.group(2),
                        "type": format_match.group(3),
                        "description": format_match.group(4),
                    }

                    if extract_formats:
                        result["formats"][format_match.group(1)] = format_field

                    # Add to fields list for compatibility
                    result["fields"].append({"key": "FORMAT", "values": format_field})

            elif line.startswith("##FILTER="):
                # process FILTER lines
                if not extract_filters and "fields" not in result:
                    continue

                # Handle different FILTER line formats
                filter_id = None
                description = ""

                # Try standard format with Description
                filter_match = re.match(
                    r"##FILTER=<ID=([^,>]+),Description=\"([^\"]+)\"",
                    line,
                )
                if filter_match:
                    filter_id = filter_match.group(1)
                    description = filter_match.group(2)
                else:
                    # Try simpler format without Description
                    filter_match = re.match(r"##FILTER=<ID=([^,>]+)>", line)
                    if filter_match:
                        filter_id = filter_match.group(1)
                        description = "No description available"
                    else:
                        # Try basic format
                        filter_match = re.match(r"##FILTER=<ID=([^>]+)>", line)
                        if filter_match:
                            filter_id = filter_match.group(1)
                            description = "No description available"

                if not filter_id:
                    logging.error("Cannot parse FILTER line: %s" % line)
                else:
                    filter_field = {
                        "ID": filter_id,
                        "description": description,
                    }

                    if extract_filters:
                        result["filters"][filter_id] = filter_field

                    # Add to fields list for compatibility
                    result["fields"].append({"key": "FILTER", "values": filter_field})

            elif line.startswith("#CHROM"):
                cols = line[1:].split()
                if extract_columns:
                    result["columns"] = cols
                # Always extract sample names (everything after FORMAT column)
                # Standard VCF columns: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [SAMPLE1] [SAMPLE2] ...
                if len(cols) > 9:  # Has samples
                    result["samples"] = cols[9:]  # Sample names start at index 9
                else:
                    result["samples"] = []
                break

    except Exception as e:
        logging.error(f"Error while reading file {filename}: {str(e)}")
        raise

    finally:
        if file_handle:
            file_handle.close()

    # Add tabix information if available
    try:
        import subprocess

        tabix_output = subprocess.check_output(
            ["tabix", "-l", filename], text=True, stderr=subprocess.PIPE
        )
        chromosomes = [
            line.strip() for line in tabix_output.split("\n") if line.strip()
        ]
        result["tabix"] = {"chromosomes": chromosomes}
    except Exception:
        result["tabix"] = None

    return result


def extract_variants(
    filename: str,
    region: Optional[str] = None,
    extract_samples: bool = True,
    sample_names: Optional[List[str]] = None,
    tabix_path: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """Extract variants from a vcf file in a specific region

    Args:
        filename: Path to VCF file (bgzipped + tabix indexed)
        region: Region to extract from, e.g. chr1:1000-2000
        extract_samples: When true, extract sample information
        sample_names: List of sample names to extract (None = all samples)
        tabix_path: Path to tabix executable

    Returns:
        List of variant records as dictionaries
    """
    if not tabix_path:
        # Try to find tabix in PATH
        tabix_path = "tabix"  # Assume it's in PATH

    command = [tabix_path, filename]

    if region:
        command.append(region)

    result = []

    p = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
    )

    header = extract_header(filename, extract_columns=True)
    columns = header.get("columns", [])

    for line in p.stdout:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        fields = line.split("\t")

        if len(fields) < 8:
            logging.error("Invalid VCF record (less than 8 fields): %s" % line)
            continue

        record = {}

        for i in range(min(len(fields), len(columns))):
            record[columns[i].lower()] = fields[i]

        # Process INFO field
        if "info" in record:
            record["info_dict"] = getInfo(record["info"])

        # Process sample data if requested
        if extract_samples and len(fields) > 9:
            record["samples"] = {}

            # Get format string
            format_keys = fields[8].split(":")

            # Extract sample data
            for i in range(9, min(len(fields), len(columns))):
                sample_name = columns[i]

                # Skip if sample_names is specified and this sample isn't in it
                if sample_names and sample_name not in sample_names:
                    continue

                sample_values = fields[i].split(":")
                sample_data = {}

                for j in range(min(len(format_keys), len(sample_values))):
                    sample_data[format_keys[j]] = sample_values[j]

                record["samples"][sample_name] = sample_data

        result.append(record)

    # Check for errors
    stderr_output = p.stderr.read()
    if stderr_output:
        logging.warning(f"Tabix stderr: {stderr_output}")

    exit_code = p.wait()
    if exit_code != 0:
        logging.error(f"Tabix exited with code {exit_code}")

    return result


def get_variant_array(
    vcf_file: str,
    region: Optional[str] = None,
    variantkeys: Optional[List[str]] = None,
    infokeys: Optional[List[str]] = None,
    samplekeys: Optional[List[str]] = None,
    tabix_path: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """Get array of values from VCF

    Args:
        vcf_file: Path to VCF file
        region: Optional region to extract from
        variantkeys: List of variant keys to extract (e.g. "chrom", "pos", "ref", "alt")
        infokeys: List of INFO fields to extract
        samplekeys: List of FORMAT fields to extract for each sample
        tabix_path: Path to tabix executable

    Returns:
        List of variant dictionaries
    """
    variants = extract_variants(
        vcf_file, region, extract_samples=samplekeys is not None, tabix_path=tabix_path
    )

    result = []

    for v in variants:
        vres = {}

        # Extract variant keys
        if variantkeys:
            for k in variantkeys:
                kl = k.lower()
                if kl in v:
                    vres[k] = v[kl]

        # Extract INFO keys
        if infokeys and "info_dict" in v:
            for k in infokeys:
                if k in v["info_dict"]:
                    vres[f"INFO_{k}"] = v["info_dict"][k]

        # Extract sample keys
        if samplekeys and "samples" in v:
            for sample, data in v["samples"].items():
                for k in samplekeys:
                    if k in data:
                        vres[f"{sample}_{k}"] = data[k]

        result.append(vres)

    return result


def writeVariantsJSON(
    vcf_file: str,
    outfile: str,
    region: Optional[str] = None,
    tabix_path: Optional[str] = None,
) -> None:
    """Extract variants and write to JSON file

    Args:
        vcf_file: Path to VCF file
        outfile: Output JSON filename
        region: Optional region to extract variants from
        tabix_path: Path to tabix executable
    """
    variants = extract_variants(vcf_file, region, tabix_path=tabix_path)

    with open(outfile, "w", encoding="utf-8") as f:
        json.dump(variants, f, indent=4)


def extractHeadersJSON(
    vcf_file: str,
    outfile: str = None,
    extract_columns: bool = True,
    extract_info: bool = True,
    extract_formats: bool = True,
    extract_filters: bool = True,
) -> Dict[str, Any]:
    """Extract VCF headers to JSON

    Args:
        vcf_file: Path to VCF file
        outfile: Optional output JSON filename. If provided, headers will be written to this file.
        extract_columns: When true, extract column headers
        extract_info: When true, extract INFO fields
        extract_formats: When true, extract FORMAT fields
        extract_filters: When true, extract FILTER fields

    Returns:
        Dictionary with header information
    """
    headers = extract_header(
        vcf_file,
        extract_columns=extract_columns,
        extract_info=extract_info,
        extract_formats=extract_formats,
        extract_filters=extract_filters,
    )

    if outfile:
        with open(outfile, "w", encoding="utf-8") as f:
            json.dump(headers, f, indent=4)

    return headers
