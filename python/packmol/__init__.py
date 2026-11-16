"""Packmol - Creates Initial Configurations for Molecular Dynamics Simulations.

The `packmol` command of this package runs the Packmol executable compiled from
Packmol.jl (https://github.com/m3g/Packmol.jl).
"""

import importlib.metadata

__version__ = importlib.metadata.version("packmol")

__all__ = ["__version__"]
