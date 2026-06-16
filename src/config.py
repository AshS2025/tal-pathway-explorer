"""
config.py
Weight profiles for pathway scoring.

A "profile" is just a pre-built list of (criterion, weight) pairs that
captures one employee/team's priorities. Calling code does:

    from config import DEFAULT_PROFILE
    scorer = WeightedPathwayScorer(DEFAULT_PROFILE)

Profiles live here (not in pathway_scoring.py) so the scoring framework
itself stays opinion-free — criteria and combinators in one place,
human preferences in another. Add new profiles by appending a new
constant below.
"""

from pathway_scoring import StepsCriterion, ThermoCriterion


# Default: balanced — same emphasis on step count as DORAnet's default
# (weight 4 on steps, weight 2 on thermo). Use this when no employee-
# specific preference is configured.
DEFAULT_PROFILE = [
    (StepsCriterion(), 4.0),
    (ThermoCriterion(), 2.0),
]
