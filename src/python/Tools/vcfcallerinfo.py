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

import itertools
import json
import logging
import os
import subprocess
import tempfile


class CallerInfo(object):
    """Class for collecting caller info and version"""

    def __init__(self):
        # callers and aligners are stored in tuples of three:
        # (caller/aligner, version, parameters)
        self.callers = []
        self.aligners = []

    def __repr__(self):
        return (
            "aligners=["
            + ",".join(["/".join(xx) for xx in self.aligners])
            + "] "
            + "callers=["
            + ",".join(["/".join(xx) for xx in self.callers])
            + "]"
        )

    def asDict(self):
        kvd = ["name", "version", "parameters"]
        return {
            "aligners": [dict(y for y in zip(kvd, x)) for x in self.aligners],
            "callers": [dict(y for y in zip(kvd, x)) for x in self.callers],
        }

    def addVCF(self, vcfname):
        """Add caller versions from a VCF
        :param vcfname: VCF file name
        """
        with tempfile.NamedTemporaryFile(delete=False) as tf:
            tf_name = tf.name

        vfh = {}
        try:
            sp = subprocess.Popen(
                "vcfhdr2json '%s' '%s'" % (vcfname, tf_name),
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            o, e = sp.communicate()

            if sp.returncode != 0:
                raise Exception(f"vcfhdr2json call failed: {o.decode('utf-8')} / {e.decode('utf-8')}")
        except Exception as ex:
            logging.error("Error running vcfhdr2json: %s", ex)
        else:
            try:
                with open(tf_name, "r", encoding="utf-8") as f:
                    vfh = json.load(f)
            except json.JSONDecodeError as ex:
                logging.error("Error parsing JSON from vcfhdr2json: %s", ex)
            except Exception as ex:
                logging.error("Error reading JSON file: %s", ex)
        finally:
            try:
                os.unlink(tf_name)
            except OSError as e:
                logging.warning(f"Error removing temp file: {e}")

        cp = ["unknown", "unknown", ""]
        gatk_callers = ["haplotypecaller", "unifiedgenotyper", "mutect"]
        sent_callers = ["haplotyper"]
        source_found = False

        for hf in vfh["fields"]:
            try:
                k = hf["key"]
                if k == "source":
                    try:
                        cp[0] = str(hf["values"])
                    except (ValueError, TypeError) as e:
                        logging.warning(f"Error parsing value: {e}")
                        continue
                        cp[0] = hf["value"]
                    if cp[0].startswith("Platypus_Version_"):
                        cp[1] = cp[0][len("Platypus_Version_") :]
                        cp[0] = "Platypus"
                    source_found = True
                elif k == "source_version":
                    try:
                        cp[1] = str(hf["values"])
                    except (ValueError, TypeError) as e:
                        logging.warning(f"Error parsing value: {e}")
                        continue
                        cp[1] = hf["value"]
                    source_found = True
                elif k == "cmdline":
                    try:
                        cp[2] = str(hf["values"])
                    except (ValueError, TypeError) as e:
                        logging.warning(f"Error parsing value: {e}")
                        continue
                        cp[2] = hf["value"]
                    source_found = True
                elif k == "platypusOptions":
                    try:
                        cp[2] = str(hf["values"])
                    except (ValueError, TypeError) as e:
                        logging.warning(f"Error parsing value: {e}")
                        continue
                        cp[2] = hf["value"]
                    source_found = True
                elif k == "octopus":
                    # octopus doesn't add a version
                    self.callers.append(["octopus", "unknown", str(hf["values"])])
                elif k.startswith("GATKCommandLine"):
                    caller = "GATK"
                    try:
                        caller += "-" + hf["values"]["ID"]
                    except (ValueError, TypeError) as e:
                        logging.warning(f"Error parsing value: {e}")
                    else:
                        version = "unknown"
                        try:
                            version = hf["values"]["Version"]
                        except (ValueError, TypeError) as e:
                            logging.warning(f"Error parsing value: {e}")
                        else:
                            options = ""
                            try:
                                options = hf["values"]["CommandLineOptions"]
                            except (ValueError, TypeError) as e:
                                logging.warning(f"Error parsing value: {e}")
                            else:
                                if any(g in caller.lower() for g in gatk_callers):
                                    self.callers.append([caller, version, options])
                elif k.startswith("SentieonCommandLine"):
                    caller = "Sentieon"
                    try:
                        caller += "-" + hf["values"]["ID"]
                    except (ValueError, TypeError) as e:
                        logging.warning(f"Error parsing value: {e}")
                    else:
                        version = "unknown"
                        try:
                            version = hf["values"]["Version"]
                        except (ValueError, TypeError) as e:
                            logging.warning(f"Error parsing value: {e}")
                        else:
                            if any(s in caller.lower() for s in sent_callers):
                                self.callers.append([caller, version])

            except (OSError, IOError) as e:
                logging.error(f"Error reading BAM file: {e}")
                continue
                pass
        if source_found:
            self.callers.append(cp)

    def addBAM(self, bamfile):
        """Extract aligner information from a BAM file
        :param bamfile: name of BAM file
        """
        sp = subprocess.Popen(
            "samtools view -H '%s'" % bamfile,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        o, e = sp.communicate()

        if sp.returncode != 0:
            raise Exception("Samtools call failed: %s / %s" % (o, e))

        for line in o.split("\n"):
            if not line.startswith("@PG"):
                continue
            try:
                # noinspection PyTypeChecker
                x = dict(y.split(":", 1) for y in line.split("\t")[1:])
            except (OSError, IOError) as e:
                logging.error(f"Error reading BAM file: {e}")
                continue
                logging.warn("Unable to parse SAM/BAM header line: %s" % line)
                continue
            cp = ["unknown", "unknown", ""]
            try:
                cp[0] = x["PN"]
            except (OSError, IOError) as e:
                logging.error(f"Error reading BAM file: {e}")
                continue
                try:
                    cp[0] = x["ID"]
                    if "-" in cp[0]:
                        cp[0] = cp[0].split("-")[0]
                except (OSError, IOError) as e:
                logging.error(f"Error reading BAM file: {e}")
                continue
                    pass
            try:
                cp[1] = x["VN"]
            except (OSError, IOError) as e:
                logging.error(f"Error reading BAM file: {e}")
                continue
                pass
            try:
                cp[2] = x["CL"]
            except (OSError, IOError) as e:
                logging.error(f"Error reading BAM file: {e}")
                continue
                pass

            self.aligners.append(cp)
