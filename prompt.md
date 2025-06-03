## 📅 Next Development Session – Road to *happy* α-release

The 2025-06-03 milestone has been delivered (see `.codex/plan_2025-06-03.md` ✔)
— automatic VCF indexing, strict typing for two key modules, and new GitHub
Actions workflows are merged.  We can now pivot towards preparing the first
α-quality build that can be exercised on the full **HG002** truth-set.

### High-level goal for the upcoming session

Ship an “alpha” artefact of the modernised codebase that:
1. Runs end-to-end on a *full-size* dataset (≈ 4 GB gzipped VCFs, whole-genome
   FASTA) without manual intervention.
2. Produces the complete deliverable set for a standard run:
   – annotated VCF **and** `.tbi`
   – summary CSV / extended CSV / JSON metrics
   – ROC TSV when `--roc` is requested
3. Passes the existing test-suite **and** a new large-dataset smoke test that
   executes only selected heavy steps in CI (behind `@pytest.mark.heavy`).

### Task list (ordered)

1. **Reference handling & template caching**
   • Confirm that `Haplo.vcfeval.runVCFEval` re-uses an SDF template across
     multiple invocations when `--scratch-prefix` is shared.  Persist the
     template in `~/.cache/happy/` to avoid repeated `rtg format` costs on
     large references.

2. **Robust vcfeval discovery**
   • Current `findVCFEval()` falls back to `rtg` on `$PATH`; add env-var
     override (`HAPPY_VCFEVAL`) and clearer error when binary is missing.

3. **Streamlined logging**
   • Replace scattered `logging.info`/`warning` calls with a central helper
     that honours `--quiet/--verbose` and timestamps.  Provide
     `--log-file <path>` CLI option.

4. **Complete type-hints phase-2**
   • Annotate `happy/qfy.py` and `Haplo/compare.py`; bump `mypy` coverage to
     these modules (`mypy.ini`, `noxfile.py`).

5. **Performance/Memory smoke-test**
   • Add `tests/heavy/test_full_dataset.py` (pytest-skipped by default) that
     references a Git-LFS placeholder or synthetic stub; executed only when
     `RUN_HEAVY=1` env-var is set.

6. **Packaging & distribution**
   • Update `pyproject.toml` metadata: classifiers, URLs.
   • Provide a minimal `Dockerfile.alpha` that pins rtg-tools and installs
     the editable package.

7. **Documentation refresh**
   • Draft `doc/alpha_release_notes.md` summarising new CLI, system
     requirements, and migration caveats.

### Stretch goals
• Remove the remaining deprecated C++/Cython shims (`src/c++` folder) from
  the default build.
• Introduce `ruff` autofix in the pre-commit chain (currently only lint).

---

_Focus on tasks 1-4 for the next session; tasks 5-7 can spill into the
following iteration once the full-dataset run succeeds._
