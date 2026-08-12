"""climt-specific exception classes."""


class ConstantNotFoundError(KeyError):
    """Raised when a required physical constant is not set in the active
    atmospheric profile.

    Subclasses KeyError so that legacy code catching KeyError still works.
    """

    def __init__(self, constant_name, profile_name=None, profile_path=None):
        self.constant_name = constant_name
        self.profile_name = profile_name
        self.profile_path = profile_path
        msg = (
            f"'{constant_name}' is not set in the current atmospheric profile. "
            f"To add it, either:\n"
            f"  1. Add it to your profile TOML under the appropriate section:\n"
            f"       {constant_name} = {{ value = ..., units = ... }}\n"
            f"  2. Set it directly: "
            f"climt.set_constant('{constant_name}', value, 'units')"
        )
        if profile_name or profile_path:
            msg += f"\nCurrent profile: {profile_name or '(custom)'} ({profile_path})"
        super().__init__(msg)

    def __str__(self):
        return self.args[0] if self.args else ""
