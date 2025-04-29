import typing as t
import logging
import sys


def get_package_logger() -> logging.Logger:
    """
    Get package logger object.
    """
    import fur_cnvkit

    return fur_cnvkit.LOGGER


def setup_logging(level: str = "INFO", name: t.Optional[str] = None) -> logging.Logger:
    """
    Configures the logging settings.
    Logs are output to stderr with a specific format and levels.
    If the same logger (by name) is passed multiple times, no duplicate handlers are added.
    """
    # 1) Determine the logging level enum from the helper
    logging_level = _get_logging_level_enum(level)

    # 2) Pick the logger: either the named one, or the package logger
    logger = logging.getLogger(name) if name else get_package_logger()

    # 3) Remove the NullHandler if present (only once)
    for handler in list(logger.handlers):
        if isinstance(handler, logging.NullHandler):
            logger.removeHandler(handler)

    # 4) Prepare the new StreamHandler + formatter
    new_handler = logging.StreamHandler(sys.stderr)
    fmt = "%(asctime)s - %(levelname)s - %(message)s"
    formatter = logging.Formatter(fmt)
    new_handler.setFormatter(formatter)

    # 5) Check for an "identical" handler already in place
    for h in logger.handlers:
        if (
            isinstance(h, logging.StreamHandler)
            and getattr(h, "stream", None) is sys.stderr
            and isinstance(h.formatter, logging.Formatter)
            and h.formatter._fmt == fmt
        ):
            # found a matching handler, so skip adding
            break
    else:
        # no existing matching handler, so safe to add
        logger.addHandler(new_handler)

    # 6) Set the logger level
    logger.setLevel(logging_level)

    return logger


def update_logger_level(logger: logging.Logger, level: str) -> None:
    """
    Update the logger level.
    """
    logging_level = _get_logging_level_enum(level)

    logger.setLevel(logging_level)
    return None


def _get_logging_level_enum(level: str) -> int:
    """
    Get the logging level enum from the logging module.
    """
    logging_level = getattr(logging, level.strip().upper())
    if not isinstance(logging_level, int):
        raise ValueError(f"Invalid log level: {level}")
    return logging_level
