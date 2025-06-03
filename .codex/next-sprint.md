    1. Green-up CI again
       • Two tests remain red after the refactor (`jeffreysCI` indexing + the CLI
         end-to-end failure triggered by a stale **site-packages** copy of
         `Tools`). Blocking α-release.
       • Decide whether we want to ship a compatibility shim (`Tools.init(*, verbose=None, **_ignore)`)
         or forcibly uninstall the old PyPI package in CI.
    2. Duplication of the *Tools* namespace
       • We now have the modernised implementation *and* whatever users may still
         have in their environment (`pip install tools`).  Consider renaming our
         package to `happy_tools` internally and expose a stub `import Tools as
         happy_tools` to avoid future collisions.
    3. Reference copy footprint
       • `~/Desktop/happy_codex/py2-hap.py` is 100 + MB; if we publish wheels or
         Docker images we should .gitignore / docker-ignore this path to prevent
         it being bundled accidentally.
    4. CLI sub-command “doctor” (environment check)
       • Planned in the remediation but not implemented yet.
    5. Heavy dataset smoke-test (@pytest.mark.heavy)
       • Needs actual dataset stubs and a GitHub Actions job behind a
         conditional.
    6. Documentation debt
       • Rewrite INSTALL section for a *non-editable* install (most users).
       • Add a short “migration guide” page that centralises the breaking changes.
    7. Performance profiling ticket
       • Now that import-time side-effects are gone, re-measure start-up time and
         memory; decide if lazy pandas import into `happy.qfy` is still needed.
    8. Future tech debt
       • Replace our internal `which()` helper with `shutil.which`.
       • Track removal of the remaining Cython shims.
    9. Release engineering
       • Bump version to `0.5.0-alpha1` once the test-suite is green again so we
         have a crisp tag for users to try.
