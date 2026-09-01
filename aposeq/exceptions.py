"""Custom exceptions for APO-SEQ."""


class ApoSeqError(Exception):
    """Base class for APO-SEQ errors."""


class ConfigError(ApoSeqError):
    """Raised when an APO-SEQ configuration is invalid."""
