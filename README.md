# sclens (omnibenchmark module)

[scLENS](https://github.com/Mathbiomed/scLENS) (Kim et al. 2024, *Nat Commun*
15:3575) as an **RDIMR** module. It takes raw counts on every FILT gene, does
its own normalization, cuts eigenvalues at the RMT edge, and keeps the signals
that survive a perturbation test.

    ./sclens.sh --output_dir out --name be1 --rawdata_h5ad be1.h5ad \
      --filtered_cellids be1_cellids.txt.gz --filtered_featureids be1_featureids.txt.gz \
      --random_seed 42

Writes `{name}_embedding.tsv` and `{name}_loadings.tsv`, plus
`{name}_sclens.json` with diagnostics (not a stage output).

- **k is derived.** `--signals robust` (the default, as published) or
  `--signals rmt` (everything above the edge, an ablation).
- **No cell is dropped.** Only the gene filter runs (`--min_cells_per_gene 15`).
- **GPU arm.** Upstream crashes on CPU when cells > genes, so `--device gpu`
  refuses to start without CUDA. Declare `requires_capabilities: [gpu]`.
- The same seed gives a byte-identical embedding.

`pca-prof: prof.sh sclens.sh` profiles under denet 0.10.3, pinned in the env,
and writes `denet-samples.jsonl` with RSS and per-PID VRAM. `prof.sh` passes
`--gpu`, which denet needs to sample VRAM. denet also discards the module's
output, so `prof.sh` saves it to `module.log` and prints it after the run.

## Measured (RTX 2000 Ada laptop, seed 42)

| dataset | cells × genes | k (rmt → robust) | wall | RSS | VRAM |
|---|---|---|---|---|---|
| duo-koh | 520 × 31364 | 34 → 34 | 2m47s | 5.1 GB | – |
| duo-zhengmix4eq | 3994 × 9460 | 14 → 12 | 2m05s | 5.1 GB | – |
| tenx-0010k | 9856 × 14632 | 74 → 57 | 9m27s | 11.6 GB | 2.1 GiB |

The first run on a host spends about 6 min installing and precompiling the
Julia packages. scLENS holds a dense cells × genes matrix, so pbmc scale is out
of reach at 64GB.
