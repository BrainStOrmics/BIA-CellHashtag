# Eval Harness

See [DESIGN.md](DESIGN.md) for architecture and [FIXTURE_SPEC.md](FIXTURE_SPEC.md)
for the JSONL schema.

## Quick start

```bash
# deterministic evals only (fast, no LLM, <30s)
pytest evals/ -m "eval and deterministic"

# full deterministic suite
pytest evals/ -m eval

# include LLM judges (requires EVAL_LLM=1 and API keys)
EVAL_LLM=1 pytest evals/ -m eval

# standalone runner (no pytest)
python evals/run_evals.py
python evals/run_evals.py --category annotation
python evals/run_evals.py --accept-baseline ann_pbmc_cd4_naive_01
```

## Regression detection

Every eval run writes `baselines/<eval_id>.run.json`. The runner compares
against `baselines/<eval_id>.baseline.json` and fails if the score drops by
more than `regression_tolerance` (per rubric).

Updating a baseline is a code-reviewed action:
```bash
python evals/run_evals.py --accept-baseline <eval_id>
git add evals/baselines/<eval_id>.baseline.json
```

## Adding a new eval

1. Add fixture(s) to `fixtures/<category>/` (≥10 records)
2. Add scorer in `scorers/<category>.py` — pure `(fixture, actual) -> float`
3. Update rubric in `rubrics/<category>.yaml` if needed
4. Generate baseline: `python evals/run_evals.py --accept-baseline <id>`
5. Document provenance in `FIXTURE_SPEC.md`
