# Labels and breakpoints

This page covers how lab MICs become benchmark labels. The code is `parse_mic` in [`data.py`](../amrbench/data.py) and `make_labels` in [`build.py`](../amrbench/build.py).

## Why log2?
Broth microdilution tests two-fold concentration steps (0.25, 0.5, 1, 2 … mg/L). On a log2 scale these steps are consecutive integers (−2, −1, 0, 1), so:
- "within ±1 dilution" is simply `|predicted − true| ≤ 1`;
- regression errors mean the same thing at low and high MICs.

## Parsing and censoring
Lab reports are often censored: the true MIC lies below the lowest or above the highest concentration tested. Every value is normalised to an integer `log2_mic` plus a censor flag:

| Lab value | `log2_mic` | `censor` | Meaning |
|---|---|---|---|
| `0.25` | −2 | `=` | exact |
| `<=0.0625` | −4 | `<=` | true ≤ 2⁻⁴ |
| `<0.125` | −4 | `<=` | true < 2⁻³, which on a dilution scale means ≤ 2⁻⁴ |
| `>64` | 6 | `>` | true > 2⁶, i.e. ≥ 2⁷ |
| `>=128` | 6 | `>` | true ≥ 2⁷ |

Unicode `≤`/`≥` are accepted.

## Breakpoints and the R/S call
The breakpoints are set in `config.breakpoints`, in mg/L. A species-specific value overrides the default. The rule is **Resistant if MIC > breakpoint**, the EUCAST "R >" convention.

With breakpoint b (in log2):

| censor | Condition | `resistant` |
|---|---|---|
| `=` | v > b | 1, otherwise 0 |
| `<=` | v ≤ b | 0 (surely S) |
| `<=` | v > b | **NaN** (unknown) |
| `>` | v ≥ b | 1 (surely R) |
| `>` | v < b | **NaN** (unknown) |

Unknown R/S labels are excluded from T5 and from ECE, but still used for MIC tasks.

> ⚠️ The pilot breakpoints are **placeholders, not EUCAST values**. Replacing them is item 3 on the [Real-data migration](Real-Data-Migration.md) checklist.

## Training target
`log2_mic_point` is a single number models can regress on: `v` for `=` and `<=`, and `v + 1` for `>`. This is a simplification. Censored-likelihood regression is an open option in the [Decision log](Decision-Log.md).

## Label file columns
`release/public/labels_trainval.tsv` and `release/private/labels_all.tsv`:

| Column | Example |
|---|---|
| `isolate_id`, `species`, `antibiotic` | `AMRB00001`, `Enterococcus faecium`, `gentamicin` |
| `mic` | `8` (original lab string) |
| `log2_mic`, `censor` | `3.0`, `=` |
| `log2_mic_point` | `3.0` |
| `breakpoint_log2` | `1.0` |
| `resistant` | `1.0` / `0.0` / empty |
| `partition` | `train` / `val` / `test` |
