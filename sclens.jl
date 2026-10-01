#!/usr/bin/env julia
"""
RDIMR module: scLENS (Kim et al. 2024, Nat Commun 15:3575).

scLENS takes raw counts and normalizes them itself (L1, log1p, z-score with
L2 cell scaling). It cuts eigenvalues at the Marchenko-Pastur/Tracy-Widom edge
fitted against a shuffled null, then keeps the signals that survive a sparse
binary perturbation test (robustness score > cos(th)). That puts it on RDIMR:
raw counts in, every FILT gene, no HVG selection. The MP fit needs the noise
genes.

Outputs
-------
{output_dir}/{name}_embedding.tsv   cell_id  PC1..PCk   (eigvec * sqrt(eigval))
{output_dir}/{name}_loadings.tsv    gene_id  PC1..PCk   (scLENS gene_basis)
{output_dir}/{name}_sclens.json     k counts, lambda, MP check, robustness (diagnostics)

k is derived. With `--signals robust` (the published method), k is the number
of signals that pass the robustness test. With `--signals rmt` (an ablation),
k is the number of eigenvalues above the TW edge. Both counts go to the JSON
either way.

GPU arm. Upstream `get_sigev` calls `cu()` with no CPU fallback when
cells > genes, so CPU runs crash there. `--device gpu` therefore refuses to
start when CUDA is not functional, rather than letting a CPU run be recorded
as the GPU arm. `--device cpu` is only for smoke tests on cells < genes.
"""

using ArgParse, CUDA, DataFrames, GZip, HDF5, Random, SparseArrays
using scLENS

function parse_args_()
    s = ArgParseSettings(description = "scLENS RDIMR module", autofix_names = false)
    @add_arg_table! s begin
        "--output_dir";          required = true
        "--name";                required = true
        "--rawdata_h5ad";        required = true; help = "rawdata H5AD (counts in layers['counts'])"
        "--filtered_cellids";    required = true; help = "gzipped list of kept cell ids"
        "--filtered_featureids"; required = true; help = "gzipped list of kept feature ids"
        "--random_seed";         required = true; arg_type = Int
        "--device";              default = "gpu"; range_tester = in(("gpu", "cpu"))
        "--signals";             default = "robust"; range_tester = in(("robust", "rmt"))
            help = "robust = signals passing the perturbation test (published); rmt = all above the TW edge (ablation)"
        "--min_cells_per_gene";  default = 15; arg_type = Int
            help = "scLENS preprocess() default. Only this gene filter runs: no cell is dropped"
        "--th";                  default = 60.0; arg_type = Float64; help = "robustness angle (degrees)"
        "--n_perturb";           default = 20; arg_type = Int
        "--p_step";              default = 0.001; arg_type = Float64
    end
    parse_args(s)
end

read_ids(path) = GZip.open(path) do io
    [strip(l) for l in eachline(io) if !isempty(strip(l))]
end

"Raw counts, cells x genes, rows/cols in the order of `cells`/`genes`."
function read_counts(path, cells, genes)
    h5open(path, "r") do f
        g = f["layers/counts"]
        enc = read_attribute(g, "encoding-type")
        enc == "csr_matrix" || error("$path: layers/counts is $enc, expected csr_matrix")
        nc, ng = read_attribute(g, "shape")
        # CSR cells x genes has the same buffers as CSC genes x cells.
        X = permutedims(SparseMatrixCSC(ng, nc, read(g, "indptr") .+ 1,
                                        read(g, "indices") .+ 1, Float32.(read(g, "data"))))
        obs, var = read(f, "obs/_index"), read(f, "var/_index")
        ci, gi = indexin(cells, obs), indexin(genes, var)
        any(isnothing, ci) && error("$(count(isnothing, ci)) filtered cell ids absent from $path")
        any(isnothing, gi) && error("$(count(isnothing, gi)) kept gene ids absent from $path")
        X = X[ci, gi]
        # Non-integer is fine (duo-koh is estimated counts): scLENS is L1 -> log1p, no count model.
        all(>=(0), nonzeros(X)) || error("$path: negative values in layers/counts; scLENS needs raw counts")
        X
    end
end

function write_tsv(path, M, ids, label)
    open(path, "w") do io
        println(io, join([label; ["PC$i" for i in 1:size(M, 2)]], '\t'))
        for (id, row) in zip(ids, eachrow(M))
            println(io, join([id; string.(row)], '\t'))
        end
    end
end

function main()
    a = parse_args_()
    println("Full command: ", join(ARGS, ' '))
    out = mkpath(a["output_dir"])
    if a["device"] == "gpu" && !CUDA.functional()
        error("--device gpu but CUDA is not functional on this host; run with --with-capability gpu on a GPU node")
    end
    Random.seed!(a["random_seed"])   # scLENS draws from the global RNG (null shuffle, perturbations)

    cells = read_ids(a["filtered_cellids"])
    genes = read_ids(a["filtered_featureids"])
    X = read_counts(a["rawdata_h5ad"], cells, genes)
    keep = vec(sum(X .!= 0, dims = 1)) .>= a["min_cells_per_gene"]
    X, genes = X[:, keep], genes[keep]
    println("  counts (cells x genes): $(size(X)), genes dropped by --min_cells_per_gene: $(count(!, keep))")
    # ponytail: fail rather than drop. A dropped cell is missing from every downstream join,
    # and scLENS's L1 step divides by the cell total.
    any(iszero, sum(X, dims = 2)) && error("cells with zero counts after the gene filter; lower --min_cells_per_gene")
    a["device"] == "cpu" && size(X, 1) > size(X, 2) &&
        error("cells > genes: upstream get_sigev needs CUDA on that path; use --device gpu")

    df = DataFrame(X, genes)
    insertcols!(df, 1, :cell => cells)
    res = scLENS.sclens(df; device_ = a["device"], th = a["th"], p_step = a["p_step"],
                        n_perturb = a["n_perturb"])
    haskey(res, :pca) || error("scLENS found no eigenvalue above the TW edge")

    n_rmt, sig = size(res[:signal_evec], 2), res[:sig_id]
    use = a["signals"] == "robust" ? sig : collect(1:n_rmt)
    isempty(use) && error("no signal passed the robustness test (n_rmt=$n_rmt); --signals rmt keeps all $n_rmt")
    emb = Matrix(res[:pca][:, 2:end])[:, use]
    load = permutedims(res[:gene_basis][use, :])
    @assert size(emb) == (length(cells), length(use)) && size(load) == (length(genes), length(use))
    all(isfinite, emb) && all(isfinite, load) || error("non-finite values in scLENS output")

    write_tsv(joinpath(out, a["name"] * "_embedding.tsv"), emb, res[:cell_id], "cell_id")
    write_tsv(joinpath(out, a["name"] * "_loadings.tsv"), load, res[:gene_id], "gene_id")
    rob = res[:robustness_scores][:rob_score]
    write(joinpath(out, a["name"] * "_sclens.json"), """
    {"n_rmt": $n_rmt, "n_robust": $(length(sig)), "k": $(length(use)), "signals": "$(a["signals"])",
     "lambda_c": $(res[:λ]), "mp_pass": $(res[:pass]), "n_cells": $(length(cells)), "n_genes": $(length(genes)),
     "device": "$(a["device"])", "random_seed": $(a["random_seed"]), "th": $(a["th"]),
     "robustness_scores": [$(join(rob, ", "))]}
    """)
    println("  n_rmt=$n_rmt n_robust=$(length(sig)) emitted k=$(length(use))")
end

main()
