import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

# Skip tests if tabix is not available
(
    pytest.skip("Skipping indexing tests: tabix not found", allow_module_level=True)
    if shutil.which("tabix") is None
    else None
)

from happy.hap import _ensure_vcf_index


@pytest.fixture
def bgzipped_vcf(tmp_path):
    # Create a minimal VCF and bgzip it
    vcf = tmp_path / "test.vcf"
    content = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    vcf.write_text(content)
    gz = tmp_path / "test.vcf.gz"
    subprocess.check_call(["bgzip", "-c", str(vcf)], stdout=open(gz, "wb"))
    return gz


def test_ensure_index_creates_tbi_and_csi(tmp_path, bgzipped_vcf):
    # Remove existing index if any
    tbi = bgzipped_vcf.with_suffix(bgzipped_vcf.suffix + ".tbi")
    csi = bgzipped_vcf.with_suffix(bgzipped_vcf.suffix + ".csi")
    if tbi.exists():
        tbi.unlink()
    if csi.exists():
        csi.unlink()
    # Create index
    _ensure_vcf_index(bgzipped_vcf, force=False)
    # Expect at least one index file
    assert tbi.exists() or csi.exists(), "Index file not created"


def test_force_recreate_index(tmp_path, bgzipped_vcf):
    tbi = bgzipped_vcf.with_suffix(bgzipped_vcf.suffix + ".tbi")
    # First create index
    _ensure_vcf_index(bgzipped_vcf)
    assert tbi.exists()
    # Write new content to .tbi to simulate old index
    tbi.write_text("dummy")
    # Force recreate
    _ensure_vcf_index(bgzipped_vcf, force=True)
    # The tbi should be overwritten (size > 5 bytes)
    assert tbi.stat().st_size > 5
