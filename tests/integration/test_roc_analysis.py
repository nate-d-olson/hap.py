import os
from pathlib import Path

import pytest

from tests.utils import get_example_dir


@pytest.mark.integration
def test_roc_analysis(tmp_path, monkeypatch):
    """Run QuantifyEngine with ROC enabled and check outputs."""
    example_dir = get_example_dir()
    truth_vcf = example_dir / "example.vcf"
    query_vcf = example_dir / "example.vcf"

    assert truth_vcf.exists(), f"Truth VCF not found: {truth_vcf}"
    assert query_vcf.exists(), f"Query VCF not found: {query_vcf}"

    # Create stub executables so dependency checks pass
    stub_dir = tmp_path / "stubs"
    stub_dir.mkdir()
    for tool in ("bgzip", "tabix"):
        stub = stub_dir / tool
        stub.write_text("#!/bin/sh\nexit 0\n")
        stub.chmod(0o755)

    monkeypatch.setenv("PATH", f"{stub_dir}{os.pathsep}{os.environ.get('PATH', '')}")

    # Import after PATH modification
    from hap_py.haplo.python_quantify import QuantifyEngine

    engine = QuantifyEngine(str(truth_vcf), str(query_vcf), enable_roc_analysis=True)

    results = engine.quantify()
    assert "roc_data" in results
    assert "bootstrap_confidence_intervals" in results

    output_prefix = tmp_path / "roc_output"
    engine.write_results(str(output_prefix))

    roc_file = output_prefix.with_suffix(".roc.tsv")
    assert roc_file.exists()
    assert roc_file.stat().st_size > 0
