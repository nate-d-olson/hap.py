Your mission for the **next coding session** is to continue the modernisation
of the `happy` repository by implementing the CLI, test, and CI tasks outlined
in `.codex/plan_2025-06-03.md`.

Focus areas (ordered):

1. **CLI** – expose the remaining vcfeval options listed in tasks C-1‒C-4 and
   ensure they are propagated down to `Haplo.compare`.
2. **VCF Indexing** – automatically bgzip & tabix the annotated VCF (task C-5).
3. **Tests** – author unit and integration coverage for the new flags, aiming
   for ≥ 90 % project coverage (tasks T-1 & T-2).
4. **Type hints** – migrate `happy/hap.py` and `Haplo/vcfeval.py` to full
   annotations and pass `mypy --strict` (see roadmap section “Type Safety”).
5. **CI migration** – draft GitHub Actions workflows replacing the legacy
   Jenkinsfile; retain parity with existing `nox` sessions.

Deliverables during the session:
• Updated source files and tests.
• Passing `nox -s tests lint type_check` locally.
• Documentation updates if any user-visible behaviour changes.

Happy hacking!  Refer back to the dated plan for acceptance criteria.
