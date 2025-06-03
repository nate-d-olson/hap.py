## 📅 Next Development Session – Road to *happy* α-release

The 2025-06-03 milestone has been delivered (see `.codex/plan_2025-06-03.md` ✔)
— automatic VCF indexing, strict typing for two key modules, and new GitHub
Actions workflows are merged.  We can now pivot towards preparing the first
α-quality build that can be exercised on the full **HG002** truth-set.

> ⚠️  **Breaking-change heads-up (2025-06-03 refactor)**

The refactor removed import-time side-effects from the *Tools* namespace and
added a `verbose` kw-arg to `Tools.init`.  Down-stream code that still calls

```python
import Tools
Tools.init()
```

continues to work, but to get INFO-level diagnostics now use

```python
Tools.init(verbose=True)
```

CI helpers and notebooks should therefore either:

1. Pass the argument explicitly *or*
2. Copy the defensive shim used in `happy.hap`:

```python
import inspect, Tools
sig = inspect.signature(Tools.init)
if "verbose" in sig.parameters:
    Tools.init(verbose=want_verbose)
else:
    Tools.init()
```

Pytest 8 also changed the behaviour of `pytest.approx`; do **not** subscript
the return object.  Instead compare directly, e.g. `assert value ==
pytest.approx(0.5)` or iterate over it.

All new documentation and examples follow these conventions.

## Task

The 2025-06-03 milestone is now complete, including:
  - Automatic VCF indexing and JSON/ROC outputs
  - Strict typing in key modules (happy.hap, Tools.init)
  - GitHub Actions workflows replacing Jenkins
  - pytest 8 compatibility fixes (removed subscripting of `pytest.approx`)

### Next Objectives (2025-06-04)
1. Prepare the first α-quality release candidate:
   - Validate end-to-end benchmarking on the full HG002 truth set (integration tests)
   - Implement and verify auto-index creation (`_ensure_vcf_index`) for all outputs (task C-5)
   - Optimize performance for large genomes and multi-threading
2. Finalize packaging and versioning:
   - Bump package version to `0.1.0a1` (alpha)
   - Update `RELEASES.md` and `setup.cfg`/`pyproject.toml` accordingly
   - Ensure `pip install .` creates the console scripts properly
3. Polish documentation and examples:
   - Update README and doc/ with HG002 usage instructions
   - Add example command-lines and expected output for full-genome benchmarking
   - Incorporate new code examples for Tools.init with `verbose`

After these tasks, update `.codex/plan_2025-06-04.md` with detailed subtasks and continue agile development toward the α-release.
