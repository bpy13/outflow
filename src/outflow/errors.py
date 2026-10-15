"""Errors raised by the outflow package."""


class NoSolutionError(RuntimeError):
    """No stationary wind solution was found for the given inputs.

    The message is written in plain language for display to users; technical details, if any,
    are in ``details``.
    """

    def __init__(self, message: str, details: str = ""):
        super().__init__(message)
        self.message = message
        self.details = details
