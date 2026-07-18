"""TE-insertion availability records.

Copy-specific TE insertion simulation is unavailable without a user-supplied TE
library, family weights, age model, and insertion preferences.
"""

from __future__ import annotations


def te_insertion_status(te_library_available: bool, insertion_preferences_available: bool) -> dict[str, object]:
    available = bool(te_library_available and insertion_preferences_available)
    return {
        "status": "available" if available else "unavailable",
        "value": None,
        "reason": (
            "User supplied both a TE library and insertion-preference model."
            if available
            else "TE insertion simulation requires a local TE library and insertion-preference model; no download is attempted."
        ),
    }
