Your mission for the **next coding session** is to finish the 2025-06-03
modernisation milestone.  The previous session completed the CLI surface and
associated tests (tasks C-1‒C-4, T-1, T-2).  The jobs that remain are:

Focus areas (ordered):

1. **VCF indexing** – automatically bgzip **and** tabix the annotated VCF that
   comes out of `Haplo.compare` (plan task C-5).  A `.vcf.gz.tbi` must always
   be present.
2. **Type hints** – migrate `happy/hap.py` and `Haplo/vcfeval.py` to full type
   annotations, enable `from __future__ import annotations`, and pass
   `mypy --strict` for these two modules (plan section “Type Safety”).
3. **CI migration** – draft GitHub Actions workflows (`python-tests.yml`,
   `lint.yml`, `type-check.yml`) that replicate the existing `nox` matrix and
   deprecate the legacy *Jenkinsfile*.
4. **Coverage gate** – ensure the new code still meets the ≥ 90 % threshold
   (task T-3) after adding indexing logic and type hints.

Deliverables
------------
• Updated source code implementing automatic indexing.
• Added or updated type annotations with `mypy --strict` clean run.
• New GitHub Actions workflow files committed under `.github/workflows/`.
• Passing `nox -s tests lint type_check` locally.
• Documentation (README, plan file) amended where user-visible behaviour
  changes.

Happy hacking!  Refer to `.codex/plan_2025-06-03.md` for exact acceptance
criteria.
