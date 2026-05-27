# Figure Generation

Plotting scripts that consume the CSVs produced by `../tests/` and emit
publication-quality PNG (300 dpi) and vector PDF files.

## Files

| Script                  | Reads                                | Outputs                                                                |
|-------------------------|--------------------------------------|------------------------------------------------------------------------|
| `plot_pareto.py`        | `../tests/data/gate_counts.csv`      | `pareto_n5/n10/n15.{png,pdf}`, `pareto_combined.{png,pdf}`             |
| `plot_compile_time.py`  | `../tests/data/compile_time.csv`     | `compile_time.{png,pdf}`                                               |
| `plot_siable.py`        | `../tests/data/siable_*.csv`         | `siable_fullrank_scaling.{png,pdf}`, `siable_low_rank_dispatch.{png,pdf}`, `siable_compile_time.{png,pdf}` |
| `_style.py`             | (imported by the three scripts)      | Shared palette, marker, and font definitions                           |

## Running

```bash
python plot_pareto.py          # state-prep Pareto figures
python plot_compile_time.py    # state-prep compile-time figure
python plot_siable.py          # SIABLE figures
```

Each script writes its output to the current directory.  All figures
are produced in both PNG and vector PDF formats.

## Visual Consistency

Every state-prep method has exactly one colour and one marker shape
across every figure, defined once in `_style.py`:

| Method                | Colour          | Marker             |
|-----------------------|-----------------|--------------------|
| low-rank              | blue            | circle             |
| PB                    | orange          | square             |
| Isometry              | green           | upward triangle    |
| low-rank+BlockZXZ     | purple          | diamond            |
| PB+BlockZXZ           | red             | downward triangle  |
| **spdmm**             | **deep pink**   | **star**           |
| **spdmm (Rust)**      | **deep pink**   | **filled plus**    |

For the SIABLE low-rank dispatch figure, the four sub-constructions
follow the same palette logic:

| SIABLE method | Colour     | Marker          | Rationale                       |
|---------------|------------|-----------------|---------------------------------|
| `state_prep`  | deep pink  | star            | uses SPDMM internally           |
| `isometry`    | green      | upward triangle | uses qiskit Isometry            |
| `spdmm_half`  | purple     | diamond         | Block-ZXZ half-isometry hybrid  |
| `full_rank`   | orange     | square          | general Block-ZXZ unitary path  |

A reader can recognise the same method by colour or marker shape on
any figure in the manuscript without consulting a separate legend.

## Recommended Manuscript Layout

- **Main text**:
  - `compile_time.pdf` (state-prep compile-time figure, responds to
    Reviewer 1 Comment 1.3 about computational overhead)
  - `pareto_combined.pdf` (single/CNOT trade-off across n = 5, 10, 15,
    responds to Reviewer 2 Comment 2.4 about single-qubit dominance)
  - `siable_fullrank_scaling.pdf` (Table III / IV evidence)

- **Supplementary appendix**:
  - `pareto_n5.pdf`, `pareto_n10.pdf`, `pareto_n15.pdf`
    (single-panel close-ups)
  - `siable_low_rank_dispatch.pdf` (Table V evidence)
  - `siable_compile_time.pdf` (SIABLE classical pre-computation cost)
