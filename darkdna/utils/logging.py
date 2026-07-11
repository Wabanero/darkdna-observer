"""Small logging facade used by CLI and library functions."""

from __future__ import annotations

import logging
import sys
from typing import Iterable


LOGGER_NAME = "darkdna"


def get_logger(name: str = LOGGER_NAME) -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
        logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    return logger


def warn_once(message: str, seen: set[str] | None = None) -> None:
    if seen is not None:
        if message in seen:
            return
        seen.add(message)
    get_logger().warning(message)


def emit_warnings(warnings: Iterable[str]) -> None:
    for warning in warnings:
        if warning:
            get_logger().warning(warning)
