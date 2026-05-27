# Tests & Benchmarks

This folder holds the correctness tests and the benchmark scripts that
produce the CSV data consumed by the plotting scripts in
`../plot_figure/`.

## Files

| File                          | Purpose                                                            |
|-------------------------------|--------------------------------------------------------------------|
| `test_siable.py`              | 8 pytest correctness tests for SIABLE (~31 s)                      |
| `_methods.py`                 | State-prep method registry + Block-ZXZ backend + benchmark helpers |
| `run_benchmark.py`            | State-prep benchmark: gate counts + compile time                   |
| `run_siable_benchmark.py`     | SIABLE benchmark: full-rank, low-rank dispatch, compile time       |
| `data/gate_counts.csv`        | Pre-computed state-prep gate counts (n = 5..15, 7 methods)         |
| `data/compile_time.csv`       | Pre-computed state-prep compile times                              |
| `data/siable_full_rank.csv`   | Pre-computed full-rank SIABLE measurements (n = 3..6)              |
| `data/siable_low_rank.csv`    | Pre-computed low-rank SIABLE dispatch table                        |

## Running the Benchmarks

### State preparation

```bash
python run_benchmark.py                    # default: n = 5..15, all methods
python run_benchmark.py --nmin 5 --nmax 10
python run_benchmark.py --methods "spdmm,Isometry,PB+BlockZXZ"
python run_benchmark.py --reps-overrides "PB:2,low-rank:2"
python run_benchmark.py --skip-timing      # only gate counts
python run_benchmark.py --fresh            # wipe existing CSVs first
```

Wall-time for the full sweep is roughly 30 minutes on a laptop.

### SIABLE

```bash
python run_siable_benchmark.py             # default: n = 3..7
python run_siable_benchmark.py --nmax 6 --reps 3
python run_siable_benchmark.py --skip-low-rank
python run_siable_benchmark.py --fresh
```

Wall-time for the default sweep is roughly 1 minute (circuit sizes
grow as O(4ⁿ) for SIABLE, so n > 7 becomes very expensive).

## Running the Correctness Tests

```bash
pytest test_siable.py -v
```

The 8 tests check the block-encoding relation, the C-NOT counts
against Table III / Table V of the paper, and the auto-dispatch logic.
