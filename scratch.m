Implementing hap.py Python 3 Migration - Getting Started
I recommend starting with this prompt and set of instructions to begin implementing the Python 3 migration for hap.py:

Initial Prompt

I need help implementing the Python 3 migration plan for hap.py, starting with fixing the build system issues for C++ components. This bioinformatics tool currently fails during the build of external C++ libraries, which is blocking further progress. Let's start by:

1. Setting up a diagnostic process to capture detailed build logs
2. Identifying critical C++ build failures from those logs
3. Creating a plan to fix those issues while working in parallel on Python-only components
4. Implementing fixes for the most critical build system issues

Please help me set up a structured approach that follows the modernization strategy in the instructions, with a focus on hap.py, qfy.py, and pre.py along with their dependencies.

Instructions to Attach
When starting this implementation session, I recommend attaching these specific instruction files:

workflow.instructions.md - For overall development workflow guidance
python3-migration.instructions.md - For Python 3 migration specifics
modernization.instructions.md - For overall modernization strategy
cython-modernization.instructions.md - For Cython-specific modernization guidance
Implementation Order
Based on the previous conversation, I suggest this implementation order:

Diagnostic phase

Capture detailed build logs
Identify specific C++ build failures
Document external library dependencies causing issues
Build system fixes

Update CMake files to address compiler issues
Fix external library integration
Create workarounds for platform-specific problems
Initial Python 3 updates (in parallel)

Start with independent modules in the Tools package
Fix bcftools.py and core utility functions
Update file handling in pre.py
Integration and testing

Test Python components with mock C++ interfaces
Verify core functionality works with the fixed build system
Add tests for the migrated components
This approach addresses both the immediate build failures while making progress on the Python 3 migration, focusing on the critical components you identified.

Would you like me to provide more specific guidance on starting with any particular component?