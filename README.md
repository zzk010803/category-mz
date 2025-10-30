# category-mz
# MW–Class Stacked Area Plotter

Tiny utility to turn a table with **m/z**, **O/C**, **H/C** into:

* class labels (Lipids / Proteins / Carbohydrates / Lignins / Condensed Aromatics / Tannins / Other),
* binned by molecular weight,
* **Origin-ready wide CSVs** (counts & %) and a stacked area PNG.

## Requirements

Python 3.8+ with `pandas`, `numpy`, `matplotlib`.

```bash
pip install pandas numpy matplotlib
```

## Quick Start

1. Edit the top of the script:

   ```python
   INPUT_PATH = r"D:\data\category-mz\NOM-0.txt"
   SAMPLE_TAG = "NOM-0"
   BIN_WIDTH_DA = 25
   ```
2. Run the script.
   Outputs are written next to your input in a folder:

   ```
   NOM-0_mw_class_outputs/
     ├─ NOM-0_mw_class_counts.csv
     ├─ NOM-0_mw_class_percent.csv
     └─ NOM-0_mw_class_stacked_area.png
   ```

## Input

Tab-delimited file with columns (case/format flexible):

* `m/z` (aka `mz`, `mass`)
* `O/C` (aka `o_c`, `oc`)
* `H/C` (aka `h_c`, `hc`)

## Classification (H/C–O/C)

* Lipids: 1.5–2.0 & 0–0.3
* Proteins: 1.5–2.2 & 0.3–0.67
* Carbohydrates: 1.5–2.2 & 0.67–1.2
* Lignins: 0.7–1.5 & 0.1–0.67
* Condensed Aromatics: 0.2–0.7 & 0–0.67
* Tannins: 0.5–1.5 & 0.67–1.2
  Else → `Other`.

## Plot in Origin (optional)

Import `*_percent.csv` → set `MW_bin` as categorical X → Area → **Stack**.

## Notes

* Change `BIN_WIDTH_DA` for different bin sizes.
* To use intensity weighting, replace `.size()` with a `.sum()` over your intensity column.

License: MIT (or adapt as needed).
