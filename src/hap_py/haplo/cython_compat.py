# Python 3 compatible imports with graceful fallback
import warnings
from typing import Any, Dict, Optional, TypeVar, Union

# Define type variables for generic return types
T = TypeVar("T")


# Define functions to import Cython modules with fallback
def import_with_fallback(
    cython_module: str, mock_attribute: str, mock_module: Optional[str] = None
) -> Any:
    """
    Import a Cython function/class with fallback to a pure Python implementation.

    Args:
        cython_module: Name of the Cython module to import
        mock_attribute: Name of the attribute to import from the mock module
        mock_module: Name of the mock module (default: 'Haplo.cython_mock')

    Returns:
        The imported object or its mock implementation
    """
    if mock_module is None:
        # Use relative import for mock module when within package
        try:
            from . import cython_mock
            mock_module = "cython_mock"
        except ImportError:
            mock_module = "hap_py.haplo.cython_mock"

    try:
        # Try to import the Cython module
        module = __import__(cython_module, fromlist=["*"])
        obj = getattr(module, mock_attribute)
        return obj
    except (ImportError, AttributeError) as e:
        # If it fails, use the mock implementation
        warnings.warn(
            f"Could not import {mock_attribute} from {cython_module}: {e}. "
            f"Using pure Python implementation from {mock_module}.",
            ImportWarning,
            stacklevel=2,
        )
        
        # Try relative import first, then absolute
        try:
            if mock_module == "cython_mock":
                from . import cython_mock as module
            else:
                module = __import__(mock_module, fromlist=["*"])
        except ImportError:
            # Fallback to local mock implementations
            from .cython import mock_cpp_internal as module
            
        obj = getattr(module, mock_attribute)
        return obj


# Import sequence utilities - try relative imports first
try:
    # Try relative import for sequence utilities within package
    complement_sequence = import_with_fallback("sequence_utils", "complement_sequence")
    reverse_complement = import_with_fallback("sequence_utils", "reverse_complement")
except ImportError:
    # Fallback to absolute imports
    complement_sequence = import_with_fallback(
        "hap_py.haplo.sequence_utils", "complement_sequence"
    )
    reverse_complement = import_with_fallback("hap_py.haplo.sequence_utils", "reverse_complement")

# Import variant processing utilities
VariantProcessor = import_with_fallback("hap_py.haplo.variant_processor", "VariantProcessor")

# Import ROC utilities
compute_roc_points = import_with_fallback("hap_py.haplo.happyroc", "compute_roc_points")

# Import chromosome sorting utilities
cmp_chromosomes = import_with_fallback("hap_py.haplo.happyroc", "cmp_chromosomes")

sort_chromosomes = import_with_fallback("hap_py.haplo.happyroc", "sort_chromosomes")


# Helper function to check if we're using Cython or pure Python
def is_using_cython() -> Dict[str, Dict[str, Union[str, bool]]]:
    """
    Check if the imported modules are Cython or pure Python.

    Returns:
        Dictionary with module information for each imported component
    """
    module_origins = {}

    for name, obj in globals().items():
        if name in [
            "complement_sequence",
            "reverse_complement",
            "VariantProcessor",
            "compute_roc_points",
            "cmp_chromosomes",
            "sort_chromosomes",
        ]:
            module = obj.__module__
            is_cython = "cython_mock" not in module
            module_origins[name] = {"module": module, "is_cython": is_cython}

    return module_origins
