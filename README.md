# sclens (omnibenchmark module)

[scLENS](https://github.com/Mathbiomed/scLENS) (Kim et al. 2024, *Nat Commun*
15:3575) as an omnibenchmark **RDIMR** module. It takes raw counts and does its
own L1 → log1p → z-score normalization with L2 cell scaling. It then cuts
eigenvalues at the Marchenko–Pastur/Tracy–Widom edge fitted against a shuffled
null, and keeps the signals that survive a sparse perturbation test. Cells come
from `filtered_cellids` and genes from `filtered_featureids`: every FILT gene,
with no HVG selection.

    ./sclens.sh --output_dir out --name be1 --rawdata_h5ad be1.h5ad \
      --filtered_cellids be1_cellids.txt.gz --filtered_featureids be1_featureids.txt.gz \
      --random_seed 42

The module writes `{name}_embedding.tsv`, `{name}_loadings.tsv` and
`{name}_sclens.json`. The JSON holds `n_rmt`, `n_robust`, λ, the MP check and the
per-signal robustness scores, and is not a declared stage output.

## Choices

- **k is derived.** `--signals robust` (the default, as published) keeps the
  signals that pass the robustness test. `--signals rmt` keeps everything above
  the TW edge, as an ablation.
- **QC.** Only scLENS's gene filter runs (`--min_cells_per_gene 15`, the
  scLENS default). No cell is dropped, so every FILT cell reaches the
  downstream joins.
- **GPU arm.** Upstream `get_sigev` calls `cu()` with no fallback when
  cells > genes. `--device gpu` (the default) refuses to start without a
  functional CUDA. `--device cpu` is for smoke tests on cells < genes only.
  The plan should declare `requires_capabilities: [gpu]`.
- Non-integer input such as duo-koh's estimated counts is accepted. scLENS has
  no count model, so only negative values are refused.
- Seeded through `Random.seed!`. On the GPU, two runs with the same seed gave
  byte-identical embeddings.

## Profiling

`pca-prof: prof.sh sclens.sh` runs the module under denet (pinned in the env,
from almost-conductor) and writes `denet-samples.jsonl` with RSS and **per-PID
VRAM**. denet ≥ 0.10 samples the GPU only when `--gpu` is passed, and
`prof.sh` passes it.

## Measured (RTX 2000 Ada laptop GPU, seed 42)

| dataset | cells × genes | n_rmt → robust | wall | RSS | VRAM |
|---|---|---|---|---|---|
| duo-koh | 520 × 31364 | 34 → 34 | 2m47s | 5.1 GB | – |
| duo-zhengmix4eq | 3994 × 9460 | 14 → 12 | 2m05s | 5.1 GB | – |
| duo-zhengmix4eq, `--min_cells_per_gene 400` (cells > genes) | 3994 × 1087 | 15 → 12 | 1m42s | 3.5 GB | 246 MiB |

Most of the wall time is Julia start-up and JIT compilation. scLENS builds a
dense cells × genes matrix and runs a full eigendecomposition of the
min(cells, genes)² Wishart matrix about 25 times, so expect it to stop scaling
well before pbmc size.

`pixi run instantiate` resolves the Julia deps (scLENS pinned by commit in
`Project.toml`, plus `Manifest.toml`). `pixi run export-env` regenerates
`envs/sclens.yml`.
