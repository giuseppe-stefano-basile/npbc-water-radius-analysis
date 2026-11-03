# Contributing

Create a branch for every methodological change.  Do not overwrite
`reference_results/`.  A contribution must include:

1. a clear statement of the estimator or convention being changed;
2. tests covering the change and a relevant limiting case;
3. output metadata recording the new convention;
4. a quantitative comparison against the reference result;
5. no machine-specific absolute paths in committed configuration files.

Before committing, run:

```bash
pytest
npbc-water check-reference reference_results
```
