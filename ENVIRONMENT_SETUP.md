# Environment Setup Documentation

## ✅ Environment Status: READY

The `happy-dev` micromamba environment has been successfully configured for hap.py development.

## Environment Details

- **Environment Name**: `happy-dev`
- **Python Version**: 3.11.12
- **Python Location**: `/Users/nolson/micromamba/envs/happy-dev/bin/python`
- **Package Installation**: Normal install in site-packages
- **Testing Framework**: pytest 8.3.5
- **RTG Tools**: Use the `rtg` executable installed in the `happy-dev` conda environment. After activating the environment, run `which rtg` to confirm its location.

## System Dependencies

Before installing hap.py in this environment, ensure that common build tools and
libraries are present. On Debian/Ubuntu systems these can be installed with:

```bash
sudo apt-get install -y build-essential python3-dev cmake zlib1g-dev libbz2-dev libboost-all-dev
```

## Required Activation Command

**ALWAYS run this before any development work:**
```bash
micromamba activate happy-dev
```

## Verification Commands

To verify the environment is working correctly:

```bash
# Check Python version and location
which python
python --version

# Verify RTG tools availability
which rtg

# Test package imports
python -c "import hap_py; print('✅ hap_py imported successfully')"
python -c "from hap_py.haplo import quantify, vcfeval; print('✅ Core modules available')"
python -c "import pytest; print(f'✅ pytest {pytest.__version__} available')"
```

## Development Workflow

1. **Activate Environment**:
   ```bash
   micromamba activate happy-dev
   ```

2. **Run Tests**:
   ```bash
   pytest tests/unit/ -v
   pytest tests/integration/ -v
   ```

3. **Code Quality**:
   ```bash
   pre-commit run --all-files
   black src/
   ruff check src/ --fix
   ```

## Project Status

Following the instructions from the prompt files:
- ✅ Environment properly activated
- ✅ Python 3.11.12 confirmed
- ✅ Package imports working
- ✅ Testing framework ready
- ✅ RTG tools available
- ✅ Ready for quantify module development

## Next Steps

The environment is ready for continuing work on the quantify module implementation as outlined in `QUANTIFY_IMPLEMENTATION_PLAN.md`.
