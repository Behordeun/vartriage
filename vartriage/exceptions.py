"""Exception and warning base classes for the vartriage package."""


class VarTriageWarning(UserWarning):
    """Base for all warnings raised by vartriage.

    Subclasses inherit from this, so users can silence everything at once:

        warnings.filterwarnings("ignore", category=VarTriageWarning)
    """
