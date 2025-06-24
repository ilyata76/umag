"""
Модуль для инструментов логирования.
"""

from .core.logger import Logger, ScopedTimer  # type: ignore # noqa


logger = Logger.instance()
