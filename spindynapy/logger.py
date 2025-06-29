"""
Logging utilities for internal performance measurement and diagnostics.

This module provides a singleton `logger` instance for use throughout the Python frontend
  as well as scoped timing tools for profiling performance-critical regions.

It wraps the C++ backend classes exposed via pybind11.

Usage:
    from spindynapy.logger import logger
    logger.add("starting simulation")
    with ScopedTimer("Simulation", always_flush=True, logger=logger):
        simulate()
"""

from .core.logger import Logger, ScopedTimer  # type: ignore # noqa


logger = Logger.instance()
"""
Global logger instance for the Python frontend.
This logger is a singleton that can be used throughout the application
"""

__all__ = ["logger", "ScopedTimer", "Logger"]
