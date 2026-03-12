#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Centralised logging configuration for the mdinterface package.

Every module in the package obtains a child logger via the standard pattern::

    import logging
    logger = logging.getLogger(__name__)

All child loggers propagate to the ``mdinterface`` parent logger, which
carries the single StreamHandler and the active level.  Call
:func:`set_verbosity` at any time to change the level for the whole package.

By default a :class:`logging.NullHandler` is installed on the parent logger
so that library log records are silently discarded when the caller has not
configured logging — following standard Python library practice.

Usage
-----
::

    import mdinterface
    mdinterface.set_verbosity("DEBUG")   # package-wide

    from mdinterface.utils.logger import set_verbosity
    set_verbosity(True)    # INFO
    set_verbosity(False)   # WARNING (silent in normal use)
    set_verbosity("DEBUG") # maximum detail
"""

import logging

_PACKAGE = "mdinterface"
_FMT = "[mdi] %(levelname)s | %(message)s"

_LEVEL_ABBREV = {
    logging.DEBUG:    "DEBG",
    logging.INFO:     "INFO",
    logging.WARNING:  "WARN",
    logging.ERROR:    "ERRR",
    logging.CRITICAL: "CRIT",
}


class _CompactFormatter(logging.Formatter):
    """Formatter that replaces level names with compact 4-char abbreviations."""

    def format(self, record):
        # Copy the record so we never mutate the original LogRecord.
        record = logging.makeLogRecord(record.__dict__)
        record.levelname = _LEVEL_ABBREV.get(record.levelno, record.levelname[:4])
        return super().format(record)


# Install a NullHandler by default so that log records from this library
# are silently discarded when the caller has not configured logging.
# This prevents the "No handlers could be found" warning and stops records
# from leaking to the root logger's lastResort handler.
logging.getLogger(_PACKAGE).addHandler(logging.NullHandler())


# Verbosity scale for small integers (0-9).
# Integers >= 10 are treated as raw logging level constants.
_VERBOSITY_SCALE = {
    0: logging.WARNING,
    1: logging.INFO,
    2: logging.DEBUG,
}


def set_verbosity(level) -> None:
    """
    Set the log level for the entire mdinterface package.

    A single StreamHandler is attached to the ``mdinterface`` root logger
    (at most once).  All child loggers (``mdinterface.build.builder``,
    ``mdinterface.io.read``, etc.) propagate to it automatically.

    Parameters
    ----------
    level : bool, int, or str
        Small integers map to a simple verbosity scale:

        - ``0`` / ``False`` -- WARNING (quiet)
        - ``1`` / ``True``  -- INFO (normal)
        - ``2``             -- DEBUG (detailed)

        Integers >= 10 are passed directly as Python logging level constants
        (e.g. ``logging.DEBUG``, ``logging.WARNING``).
        Strings are resolved by name (``"DEBUG"``, ``"INFO"``, ``"WARNING"`` ...).
    """
    if isinstance(level, bool):
        level = logging.INFO if level else logging.WARNING
    elif isinstance(level, int) and level < 10:
        level = _VERBOSITY_SCALE.get(level, logging.DEBUG)
    elif isinstance(level, str):
        level = getattr(logging, level.upper(), logging.INFO)

    parent = logging.getLogger(_PACKAGE)

    # Add a StreamHandler exactly once (NullHandler is not a StreamHandler).
    if not any(isinstance(h, logging.StreamHandler) for h in parent.handlers):
        handler = logging.StreamHandler()
        handler.setFormatter(_CompactFormatter(_FMT))
        parent.addHandler(handler)

    parent.setLevel(level)
    parent.propagate = False  # don't bubble up to the root logger

    # Make sure all existing child loggers propagate to the parent
    # and do not filter records themselves.
    for name, obj in logging.Logger.manager.loggerDict.items():
        if name.startswith(_PACKAGE + ".") and isinstance(obj, logging.Logger):
            obj.handlers.clear()
            obj.setLevel(logging.NOTSET)
            obj.propagate = True


_HEADER_WIDTH = 44  # total character width of header and subheader lines


def log_header(log: logging.Logger, title: str,
               level: int = logging.INFO) -> None:
    """
    Emit a fixed-width section header via *log*, preceded by a blank line.

    The title is left-anchored after ``===`` and ``=`` fills the rest to
    :data:`_HEADER_WIDTH` characters, so headers align regardless of title
    length.

    Example
    -------
    ::

        from mdinterface.utils.logger import log_header
        log_header(logger, "Build")
        # [mdi] INFO |
        # [mdi] INFO | ===  Build  ================================
    """
    log.log(level, "")
    left = f"=== {title} "
    fill = "=" * (_HEADER_WIDTH - len(left))
    log.log(level, "%s%s", left, fill)


def log_banner(log: logging.Logger, *lines: str,
               level: int = logging.INFO) -> None:
    """
    Emit a prominent multi-line banner with ``=`` borders.

    Each line is centred within :data:`_HEADER_WIDTH` characters.  Intended
    for the top-level startup message of a major component (e.g. SimCell).

    Example
    -------
    ::

        log_banner(logger, "mdinterface :: SimCell", "version 1.5.0")
        # [mdi] INFO | ============================================
        # [mdi] INFO |        mdinterface :: SimCell
        # [mdi] INFO |              version 1.5.0
        # [mdi] INFO | ============================================
    """
    bar = "=" * _HEADER_WIDTH
    log.log(level, bar)
    for line in lines:
        log.log(level, line.center(_HEADER_WIDTH))
    log.log(level, bar)


def log_subheader(log: logging.Logger, title: str,
                  level: int = logging.INFO) -> None:
    """
    Emit a fixed-width sub-section header via *log*.

    Same width and left-anchor style as :func:`log_header` but uses ``-``
    as the fill character to indicate a lower level of hierarchy.

    Example
    -------
    ::

        log_subheader(logger, "Layer [1/4]")
        # [mdi] INFO | --  Layer [1/4]  ----------------------------
    """
    log.log(level, "")
    left = f"-- {title} "
    fill = "-" * max(2, _HEADER_WIDTH - len(left))
    log.log(level, "%s%s", left, fill)
