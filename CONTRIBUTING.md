# Contributing to vartriage

Thanks for considering a contribution. Here's how to get set up, run tests, and open a PR.

## Development Setup

Clone the repository and install in editable mode with development dependencies:

```bash
git clone https://github.com/Behordeun/vartriage.git
cd vartriage
pip install -e ".[dev]"
```

This installs pytest, hypothesis, pytest-cov, and mypy alongside the core package.

### Accelerated Backends

For the optional fast backends (polars, pyranges) and PDF support:

```bash
pip install -e ".[all]"
```

This pulls in polars, pyranges, and reportlab.

## Running Tests

Run the full test suite with:

```bash
pytest
```

A passing run looks like:

```text
tests/ ... 383 passed in Xs
```

The project uses [Hypothesis](https://hypothesis.readthedocs.io/) for property-based testing with three profiles:

| Profile   | Max Examples | Usage                          |
| --------- | ------------ | ------------------------------ |
| `dev`     | 50           | Default for local development  |
| `ci`      | 500          | Used in CI pipelines           |
| `debug`   | 10           | Quick iteration when debugging |

Switch profiles via the `HYPOTHESIS_PROFILE` environment variable:

```bash
HYPOTHESIS_PROFILE=ci pytest
```

To skip slow performance benchmarks during local development:

```bash
pytest -m "not slow"
```

## Type Checking

Strict mypy is enforced:

```bash
mypy --strict
```

This uses the config in `pyproject.toml` (Python 3.10 target, strict mode). A clean run means zero errors.

## Code Style

We use **Black** for formatting and **ruff** for linting. Before opening a PR:

```bash
black .
ruff check .
```

Fix any auto-fixable lint issues with:

```bash
ruff check --fix .
```

## Branch Naming

Use the following prefixes for your branches:

- `feature/` — new functionality (e.g., `feature/streaming-reports`)
- `fix/` — bug fixes (e.g., `fix/score-loader-nan-handling`)
- `docs/` — documentation changes (e.g., `docs/update-api-reference`)

## Pull Request Process

1. Create a branch from `main` using the naming convention above.
2. Make your changes, keeping commits focused on a single logical change.
3. Make sure CI passes:
   - `pytest` (full suite)
   - `mypy --strict` (zero errors)
   - Black + ruff (formatting/lint)
4. Open a pull request against `main`.
5. PRs require at least one approval before merging.
6. Keep the PR description short: what changed, what you tested, open questions if any.

## Project Structure

```text
vartriage/
├── annotation/       # Variant annotation engine
├── classification/   # ACMG classification and combining rules
├── models/           # Data models (Variant, AnnotatedVariant, etc.)
├── prioritization/   # Scoring and prioritization engine
├── reporting/        # Report generation (JSON, CSV, PDF)
├── _internal/        # Internal utilities
├── protocols.py      # Protocol interfaces
├── exceptions.py     # Warning and exception hierarchy
├── cli.py            # Command-line interface
└── py.typed          # PEP 561 marker
```
