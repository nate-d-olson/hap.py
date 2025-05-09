hap.py is a bioinformatics tool for benchmarking small variant calls.
The tool is no longer maintained with outdated dependencies and used 
python2, which is no longer supported.
The tool is difficult to install limiting usability.
Despite these limitations it is a popular bioinformatics tool.

This is a fork of the developer repository.

Currently in the process of converting the tool to use python3 and update depdencies to improve usability and facilitate future development. After updates the codebase should follow python best practices and modern python package structure and instalation. Black should be used for linting.

Currently the install process is handled by the `install.py` script. Which accepts a build directory as a command line argument. Use a temp directory when building and testing. Do not make edits to the build directory as they will be lost next time the install script is run. The tool uses cython. CMake is used for installing c++ dependencies and libraries and building source code. After the install completes the install script tests the install using a series of bash and python scripts in the `src/sh` directory which are executed by the `src/sh/run_tests.sh` script. The test scripts should be run without failing on error so that all errors can be indentified. As the script contain multiple tests with ambiguous names the standard error and standard out should be recorded and used to determine the individual failing tests and specific part of the codebase that is failing the test.

Plan for updating tool.

- Convert code based from python2 to python3. Python 2to3 package has been used for automated convertion. Syntax errors and other issues remain, primarily related to changes to function arguments for python dependencies.
- Run test and document errors. After collecting information on all the failing tests, review codebase and begin debugging common syntax errors. Interatively test and debug the code until all the tests pass.
- Restructure the codebase to follow modern python base practices for python packages, syntax, and installation processes.
- Make code more robust by adding typehints and assertions as appropriate.
- Convert current tests to a formal testing framework, e.g. pytest, unittest, or whatever is deemed more approriate after code review.
- Update c++ and external libraries as appropriate.
- Clean-up codebase by removing old code and functionality that is no used and artifact of initial tool development and update process.