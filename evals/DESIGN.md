# Eval Layer — CellHashtag (OpenAI Harness Standard)

This directory is the eval harness for CellHashtag. It follows the OpenAI Harness philosophy:
**standardized fixtures, rubric-driven scoring, regression detection, and strict non-interference
with production code.** Nothing under `evals/` is imported by `src/cellhashtag/`. Ever.

## Directory Layout

```
evals/
├── DESIGN.md                 ← you are here
├── FIXTURE_SPEC.md           ← JSONL schema + per-category examples
├── conftest.py               ← pytest markers, fixture loaders, skip rules
├── run_evals.py              ← standalone runner (no pytest required)
│
├── fixtures/                 ← JSONL input/expected pairs, checked into git
│   ├── cluster_quality/      ← (adata_pickle_path, expected_asw_range, expected_k)
│   ├── annotation/           ← (markers, tissue, gold_cell_type)
│   ├── mcts/                 ← (root_state, expected_best_action, expected_depth)
│   ├── tool_use/             ← (query, expected_tool, expected_args)
│   └── ca_aa/                ← (annotation_state, expected_next_action)
│
├── rubrics/                  ← YAML rubrics per category (dims, scale, thresholds)
│   ├── cluster_quality.yaml
│   ├── annotation.yaml
│   ├── mcts.yaml
│   ├── tool_use.yaml
│   └── ca_aa.yaml
│
├── scorers/                  ← deterministic scorers (pure functions, no LLM)
│   ├── cluster_quality.py
│   ├── annotation.py
│   ├── mcts.py
│   ├── tool_use.py
│   └── ca_aa.py
│
├── judges/                   ← LLM-judged scorers (adversarial, calibrated)
│   └── adversarial_annotation.py
│
└── baselines/                ← snapshot JSON per eval, used for regression detection
    └── <eval_id>.baseline.json
```

## 1. Fixture Format

See `FIXTURE_SPEC.md`. Summary: JSONL, one record per line, schema:

```json
{"id": "ann_pbmc_cd4_t_01", "category": "annotation",
 "input": {...}, "expected": {...}, "metadata": {...}}
```

Fixtures are **immutable reference data**. Gold labels come from:
- Published PBMC datasets (10x Genomics, Zheng et al. 2017)
- CellMarker-curated canonical marker lists
- Manually reviewed annotations by a domain expert (for small calibration sets)

## 2. Rubric Design

Each rubric YAML has: `dimensions`, `scale` (0-1 or discrete), `weights`,
`pass_threshold`, `fail_threshold`, `regression_tolerance`.
See `rubrics/*.yaml` for the five category rubrics.

## 3. Integration

| Trigger | Mechanism | What runs |
|---|---|---|
| Local dev | `pytest evals/ -m eval` | all deterministic evals |
| Local dev, fast | `pytest evals/ -m "eval and not llm"` | no-LLM subset (<30s) |
| CI (PR) | `python evals/run_evals.py --ci` | deterministic + small LLM judge set |
| Nightly | `python evals/run_evals.py --full` | full LLM judge + regression check |
| Claude hook | `.claude/hooks/post-edit.sh` | fast subset on every `src/cellhashtag/` edit |

Pytest markers (defined in `conftest.py`):
- `eval` — all evals
- `deterministic` — pure-function, no LLM
- `llm` — requires LLM call (skipped without `EVAL_LLM=1`)
- `slow` — MCTS convergence / full CA-AA loop (>60s)
- `regression` — compares to `baselines/` snapshot

## 4. Deterministic vs LLM-Judged Evals

**Deterministic** (preferred whenever possible):
- Cluster silhouette, modularity, cluster count
- Annotation exact match / Jaccard vs gold set
- MCTS best-action match vs gold policy
- Tool-use argument exact match
- CA-AA state transition match vs gold policy

**LLM-judged** (only when no deterministic signal exists):
- Annotation free-text reasoning quality
- Adversarial cross-validator score (dim 3 in reward function)

LLM judges are **always adversarial + calibrated**:
- Judge uses a *different* model/temperature than the system under test
- Judge prompt is adversarial ("find why this is wrong")
- Judge output is calibrated against a 50-record human-labeled calibration set
  (`fixtures/annotation/calibration.jsonl`); judge is rejected if agreement < 0.75 κ
- Judge score is combined with deterministic dims (never used alone)
  — see rubrics/annotation.yaml weight split (≥60% deterministic)

## 5. Regression Detection

Each eval run writes a snapshot to `baselines/<eval_id>.run.json`. The runner
compares against `baselines/<eval_id>.baseline.json`:

- **Pass**: score ≥ pass_threshold
- **Regression**: score dropped by > regression_tolerance from baseline
- **Drift**: score within tolerance but directionally worse for 3 consecutive runs
  (tracked in `baselines/<eval_id>.history.jsonl`)

Regressions **fail the CI job** and require manual baseline update via:
```
python evals/run_evals.py --accept-baseline <eval_id>
```

Baseline updates are code-reviewed like production changes.

## 6. Non-Interference Guarantee

- `evals/` is a sibling of `src/`, never on `sys.path` of production code
- Production code has **zero** `evals` imports (enforced by `conftest.py`
  import-linter check)
- Fixtures reference production via `src.cellhashtag.*` imports inside
  eval test functions only
- `pyproject.toml` excludes `evals/` from package wheels

## 7. Adding a New Eval

1. Add fixture(s) under `fixtures/<category>/` (≥10 records for statistical power)
2. Add scorer in `scorers/<category>.py` — pure function `(fixture, actual) -> Score`
3. If rubric changes, edit `rubrics/<category>.yaml`
4. Add pytest test in `tests/test_<category>.py` with `@pytest.mark.eval`
5. Generate baseline: `python evals/run_evals.py --accept-baseline <eval_id>`
6. Document fixture provenance in `FIXTURE_SPEC.md`

## 8. Anti-Gaming / Anti-Sycophancy in Evals

The eval layer itself applies the same anti-sycophancy principles as the
production adversarial evaluation:

- **Hard evidence first**: ≥60% of every rubric weight is deterministic
- **Independent judges**: LLM judge uses different temperature than system under test
- **Adversarial framing**: judge prompts ask "why is this wrong", not "how good is this"
- **Calibration set**: judge agreement with human labels is measured; low-κ judges
  are rejected, not deployed
- **Variance penalty**: near-tie top-2 predictions are penalized (see rubrics/mcts.yaml)
- **Overconfidence penalty**: predicted confidence vs empirical accuracy is tracked;
  ECE (Expected Calibration Error) > 0.15 fails the eval
