#!/usr/bin/env python3
"""
Test script for checking mock implementations of Cython modules

This script tests if the mock implementations for Cython modules
work correctly when the actual C++ components are not available.

Usage:
    python3 test_cython_mock.py
"""

import os
import sys
import logging
import importlib

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

def test_mock_implementations():
    """Test the mock implementations of Cython modules"""
    
    # Force using mock implementation
    os.environ["HAPLO_USE_MOCK"] = "1"
    
    try:
        # Add the src directory to sys.path
        src_python_dir = os.path.join(os.getcwd(), 'src', 'python')
        if os.path.exists(src_python_dir):
            if src_python_dir not in sys.path:
                sys.path.insert(0, src_python_dir)
                logging.info(f"Added {src_python_dir} to Python path")
        else:
            logging.error(f"Directory not found: {src_python_dir}")
            
        logging.info(f"Python path: {sys.path}")
        
        # Try importing the module with the mock implementation
        import Haplo.cython
        from Haplo.cython import mock__internal, mock_cpp_internal
        from Haplo.cython import USING_MOCK
        
        if not USING_MOCK:
            logging.error("USING_MOCK flag is False, but should be True")
            return False
            
        # Test basic functionality
        test_result = mock__internal.test_string_handling()
        logging.info(f"Mock _internal.test_string_handling returned: {test_result}")
        
        test_result = mock_cpp_internal.test_string_handling()
        logging.info(f"Mock cpp_internal.test_string_handling returned: {test_result}")
        
        logging.info("Mock implementation working correctly!")
        return True
        
    except ImportError as e:
        logging.error(f"Failed to import mock modules: {e}")
        return False
    except Exception as e:
        logging.error(f"Error testing mock implementation: {e}")
        return False

if __name__ == "__main__":
    success = test_mock_implementations()
    sys.exit(0 if success else 1)
