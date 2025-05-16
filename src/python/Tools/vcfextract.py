#!/usr/bin/env python33
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

import os
import logging
import subprocess
import tempfile
import gzip
import re
import time
import json
from typing import Dict, List, Union, Any, Optional, Tuple


def field(val: str) -> Union[int, float, str, List[Any]]:
    """ extract field into result, guess type """
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


def extract_header(filename: str, extract_columns: bool = True, 
                  extract_info: bool = False, extract_formats: bool = False) -> Dict[str, Any]:
    """Extract header information from a VCF file

    Args:
        filename: Path to VCF file
        extract_columns: When true, extract column headers
        extract_info: When true, extract INFO fields
        extract_formats: When true, extract FORMAT fields

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

    fh = None
    try:
        if filename.endswith(".gz"):
            fh = gzip.open(filename, "rt", encoding='utf-8')
        else:
            fh = open(filename, "rt", encoding='utf-8')

        for l in fh:
            l = l.strip()
            if l.startswith("##INFO="):
                # process INFO lines
                if not extract_info:
                    continue
                
                m = re.match(r"##INFO=<ID=([^,>]+),Number=([^,>]+),Type=([^,>]+),Description=\"([^\"]+)\"", l)
                if not m:
                    logging.error("Cannot parse INFO line: %s" % l)
                else:
                    result["info"][m.group(1)] = {"id": m.group(1),
                                                "number": m.group(2),
                                                "type": m.group(3),
                                                "description": m.group(4)}
            elif l.startswith("##FORMAT="):
                # process FORMAT lines
                if not extract_formats:
                    continue

                m = re.match(r"##FORMAT=<ID=([^,>]+),Number=([^,>]+),Type=([^,>]+),Description=\"([^\"]+)\"", l)
                if not m:
                    logging.error("Cannot parse FORMAT line: %s" % l)
                else:
                    result["formats"][m.group(1)] = {"id": m.group(1),
                                                  "number": m.group(2),
                                                  "type": m.group(3),
                                                  "description": m.group(4)}
            elif l.startswith("#CHROM"):
                if extract_columns:
                    cols = l[1:].split()
                    result["columns"] = cols
                break

    except Exception as e:
        logging.error("Error while reading file %s: %s" % (filename, str(e)))
        raise

    finally:
        if fh:
            fh.close()

    return result


def extract_variants(filename: str, 
                    region: Optional[str] = None, 
                    extract_samples: bool = True, 
                    sample_names: Optional[List[str]] = None,
                    tabix_path: Optional[str] = None) -> List[Dict[str, Any]]:
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

    p = subprocess.Popen(command,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True)
    
    header = extract_header(filename, extract_columns=True)
    columns = header.get("columns", [])

    for l in p.stdout:
        l = l.strip()
        if not l or l.startswith("#"):
            continue

        fields = l.split("\t")
        
        if len(fields) < 8:
            logging.error("Invalid VCF record (less than 8 fields): %s" % l)
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


def get_variant_array(vcf_file: str, 
                     region: Optional[str] = None, 
                     variantkeys: Optional[List[str]] = None, 
                     infokeys: Optional[List[str]] = None,
                     samplekeys: Optional[List[str]] = None,
                     tabix_path: Optional[str] = None) -> List[Dict[str, Any]]:
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
    variants = extract_variants(vcf_file, region, extract_samples=samplekeys is not None, tabix_path=tabix_path)
    
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


def writeVariantsJSON(vcf_file: str, 
                     outfile: str, 
                     region: Optional[str] = None,
                     tabix_path: Optional[str] = None) -> None:
    """Extract variants and write to JSON file
    
    Args:
        vcf_file: Path to VCF file
        outfile: Output JSON filename
        region: Optional region to extract variants from
        tabix_path: Path to tabix executable
    """
    variants = extract_variants(vcf_file, region, tabix_path=tabix_path)
    
    with open(outfile, "w", encoding='utf-8') as f:
        json.dump(variants, f, indent=4)
