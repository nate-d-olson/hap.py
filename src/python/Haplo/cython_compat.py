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
        mock_module = "Haplo.cython_mock"

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
        module = __import__(mock_module, fromlist=["*"])
        obj = getattr(module, mock_attribute)
        return obj


# Import sequence utilities
complement_sequence = import_with_fallback(
    "Haplo.sequence_utils", "complement_sequence"
)

reverse_complement = import_with_fallback("Haplo.sequence_utils", "reverse_complement")

# Import variant processing utilities
VariantProcessor = import_with_fallback("Haplo.variant_processor", "VariantProcessor")

# Import ROC utilities
compute_roc_points = import_with_fallback("Haplo.happyroc", "compute_roc_points")

# Import chromosome sorting utilities
cmp_chromosomes = import_with_fallback("Haplo.happyroc", "cmp_chromosomes")

sort_chromosomes = import_with_fallback("Haplo.happyroc", "sort_chromosomes")


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
