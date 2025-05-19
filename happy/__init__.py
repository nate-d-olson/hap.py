# Top-level happy package to expose src/python/happy for development imports
import os
# Extend package path to include source directory
__path__.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src', 'python', 'happy'))
)