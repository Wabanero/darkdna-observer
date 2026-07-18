# Phase 2 and Phase 3 verification

Verification completed on 2026-07-18 with the project `genetichyper`
environment.

## Automated tests

- `pytest -m "not integration" -q`: 62 passed.
- `pytest tests/test_integration_toy_pipeline.py -q`: 1 passed in 178.7 s.
- `pytest tests/test_integration_yeast_pipeline.py -q`: 1 passed in 30.4 s
  when rerun together with the primitive-labeler regression tests.
- `pytest tests/test_integration_arabidopsis_pipeline.py -q`: 1 passed in
  157.4 s.
- `python -m compileall -q darkdna tests`: passed.
- `git diff --check`: passed.

The integration tests execute the complete legacy Mode A CLI path from window
generation through reports and browser tracks. Dedicated unit and smoke tests
cover the severe null registry and calibration, evolutionary null simulator,
native-versus-randomized benchmark, quantity/copy/presence/length/spacing
features, sequence-indifference perturbations, Mode B scoring, and Mode A-B
comparison.

## Scientific promotion rule

Primitive labels remain candidate hypotheses. A label is eligible for candidate
promotion only when at least three named null models are block-calibrated, the
configured agreement threshold is met, no null-family conflict is present, and
the conservative matched-null z-score passes its threshold. Missing or
conflicting null evidence is represented explicitly and retains the result as a
screening-only hypothesis.

Mode B perturbations are model-based sensitivity experiments. They do not, on
their own, establish biological causality, function, or selected-effect status.
