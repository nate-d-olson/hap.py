#
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/sequencing/licenses/blob/master/Simplified-BSD-License.txt

import contextlib
import json
import logging
import os
import subprocess
import tempfile


class CallerInfo:
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
        tf = tempfile.NamedTemporaryFile(delete=False)
        tf.close()
        vfh = {}
        try:
            sp = subprocess.Popen(
                f"vcfhdr2json '{vcfname}' '{tf.name}'",
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            o, e = sp.communicate()

            if sp.returncode != 0:
                raise Exception(f"vcfhdr2json call failed: {o} / {e}")

            vfh = json.load(open(tf.name))
        finally:
            with contextlib.suppress(Exception):
                os.unlink(tf.name)

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
                    except Exception:
                        cp[0] = hf["value"]
                    if cp[0].startswith("Platypus_Version_"):
                        cp[1] = cp[0][len("Platypus_Version_") :]
                        cp[0] = "Platypus"
                    source_found = True
                elif k == "source_version":
                    try:
                        cp[1] = str(hf["values"])
                    except Exception:
                        cp[1] = hf["value"]
                    source_found = True
                elif k == "cmdline":
                    try:
                        cp[2] = str(hf["values"])
                    except Exception:
                        cp[2] = hf["value"]
                    source_found = True
                elif k == "platypusOptions":
                    try:
                        cp[2] = str(hf["values"])
                    except Exception:
                        cp[2] = hf["value"]
                    source_found = True
                elif k == "octopus":
                    # octopus doesn't add a version
                    self.callers.append(["octopus", "unknown", str(hf["values"])])
                elif k.startswith("GATKCommandLine"):
                    caller = "GATK"
                    with contextlib.suppress(Exception):
                        caller += "-" + hf["values"]["ID"]

                    version = "unknown"
                    with contextlib.suppress(Exception):
                        version = hf["values"]["Version"]

                    options = ""
                    with contextlib.suppress(Exception):
                        options = hf["values"]["CommandLineOptions"]

                    if any(g in caller.lower() for g in gatk_callers):
                        self.callers.append([caller, version, options])
                elif k.startswith("SentieonCommandLine"):
                    caller = "Sentieon"
                    with contextlib.suppress(Exception):
                        caller += "-" + hf["values"]["ID"]

                    version = "unknown"
                    with contextlib.suppress(Exception):
                        version = hf["values"]["Version"]

                    options = ""
                    if any(s in caller.lower() for s in sent_callers):
                        self.callers.append([caller, version])

            except Exception:
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
            raise Exception(f"Samtools call failed: {o} / {e}")

        for line in o.split("\n"):
            if not line.startswith("@PG"):
                continue
            try:
                # noinspection PyTypeChecker
                x = dict(y.split(":", 1) for y in line.split("\t")[1:])
            except Exception:
                logging.warn("Unable to parse SAM/BAM header line: %s" % line)
                continue
            cp = ["unknown", "unknown", ""]
            try:
                cp[0] = x["PN"]
            except Exception:
                try:
                    cp[0] = x["ID"]
                    if "-" in cp[0]:
                        cp[0] = cp[0].split("-")[0]
                except Exception:
                    pass
            with contextlib.suppress(Exception):
                cp[1] = x["VN"]

            with contextlib.suppress(Exception):
                cp[2] = x["CL"]

            self.aligners.append(cp)
