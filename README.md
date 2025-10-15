# Elliptic Curve Discrete Logarithm Experiments

This repository contains scripts and datasets for studying the performance of
classical algorithms that solve the elliptic-curve discrete logarithm problem
(ECDLP) on small prime fields. Two entry points are provided:

- `ecdlp_experiments.py` generates fresh trial data for brute force,
  baby-step/giant-step, Pollard's rho, and Pollard's kangaroo solvers across the
  six sample curves.
- `ecdlp_analysis.py` aggregates previously recorded trials and produces a
  Markdown summary with descriptive statistics and (optionally) plots.

## Downloading the data and scripts

You can retrieve the files in several ways:

1. **Clone the repository** (recommended if you plan to rerun experiments):
   ```bash
   git clone <repository-url>
   cd <repository-directory>
   ```

2. **Download a ZIP archive** directly from the hosting platform (e.g. GitHub):
   - Navigate to the repository page in your browser.
   - Choose **Code → Download ZIP**.
   - Extract the archive locally to access the CSVs, scripts, and report.

3. **Fetch individual files** via your browser:
   - Open the file in the repository view.
   - Click the **Raw** button and use your browser's *Save As* option to store
     the file locally.

After downloading, you can run the experiments or analysis with:

```bash
python ecdlp_experiments.py
python ecdlp_analysis.py
```

Both scripts are dependency-free; they rely only on the Python standard library.
