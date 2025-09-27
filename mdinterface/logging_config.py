"""Logging configuration for mdinterface package."""

import logging
import logging.config
import sys
from pathlib import Path
from typing import Optional


def setup_logging(
    level: str = "INFO",
    log_file: Optional[Path] = None,
    format_string: Optional[str] = None,
) -> None:
    """Configure logging for mdinterface package.

    Parameters
    ----------
    level : str, default="INFO"
        Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    log_file : Path, optional
        Path to log file. If None, logs only to console
    format_string : str, optional
        Custom format string for log messages

    Examples
    --------
    >>> setup_logging(level="DEBUG")
    >>> setup_logging(level="INFO", log_file=Path("mdinterface.log"))
    """
    if format_string is None:
        format_string = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

    # Configure root logger
    logging_config = {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "standard": {
                "format": format_string,
                "datefmt": "%Y-%m-%d %H:%M:%S",
            },
            "detailed": {
                "format": (
                    "%(asctime)s - %(name)s - %(levelname)s - "
                    "%(filename)s:%(lineno)d - %(funcName)s - %(message)s"
                ),
                "datefmt": "%Y-%m-%d %H:%M:%S",
            },
        },
        "handlers": {
            "console": {
                "class": "logging.StreamHandler",
                "level": level,
                "formatter": "standard",
                "stream": sys.stdout,
            },
        },
        "loggers": {
            "mdinterface": {
                "level": level,
                "handlers": ["console"],
                "propagate": False,
            },
        },
        "root": {
            "level": "WARNING",
            "handlers": ["console"],
        },
    }

    # Add file handler if log_file is specified
    if log_file is not None:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)

        logging_config["handlers"]["file"] = {
            "class": "logging.handlers.RotatingFileHandler",
            "level": level,
            "formatter": "detailed",
            "filename": str(log_file),
            "maxBytes": 10_000_000,  # 10MB
            "backupCount": 5,
        }
        logging_config["loggers"]["mdinterface"]["handlers"].append("file")

    logging.config.dictConfig(logging_config)


def get_logger(name: str) -> logging.Logger:
    """Get a logger instance for mdinterface modules.

    Parameters
    ----------
    name : str
        Name of the logger (usually __name__)

    Returns
    -------
    logging.Logger
        Configured logger instance

    Examples
    --------
    >>> logger = get_logger(__name__)
    >>> logger.info("This is an info message")
    """
    return logging.getLogger(name)


class ProgressLogger:
    """Logger for tracking progress of long-running operations."""

    def __init__(self, name: str, total_steps: int):
        """Initialize progress logger.

        Parameters
        ----------
        name : str
            Name of the operation being tracked
        total_steps : int
            Total number of steps in the operation
        """
        self.logger = get_logger(f"mdinterface.progress.{name}")
        self.name = name
        self.total_steps = total_steps
        self.current_step = 0

    def step(self, message: str = "") -> None:
        """Log a step in the progress.

        Parameters
        ----------
        message : str, optional
            Additional message to include with the step
        """
        self.current_step += 1
        progress_pct = (self.current_step / self.total_steps) * 100

        log_message = f"{self.name}: Step {self.current_step}/{self.total_steps} ({progress_pct:.1f}%)"
        if message:
            log_message += f" - {message}"

        self.logger.info(log_message)

    def complete(self, message: str = "Completed successfully") -> None:
        """Log completion of the operation.

        Parameters
        ----------
        message : str, default="Completed successfully"
            Completion message
        """
        self.logger.info(f"{self.name}: {message}")


# Initialize default logging on import
setup_logging()
