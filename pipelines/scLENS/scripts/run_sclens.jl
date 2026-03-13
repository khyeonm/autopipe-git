import Pkg
Pkg.activate("/pipeline")

using scLENS
using CUDA: has_cuda
using CSV
using DataFrames
using CairoMakie

# Arguments: input_csv pca_out umap_out h5ad_out plot_out
args = ARGS
input_csv = args[1]
pca_out   = args[2]
umap_out  = args[3]
h5ad_out  = args[4]
plot_out  = args[5]

println("Loading data from: $input_csv")
ndf = scLENS.read_file(input_csv)

println("Running QC / preprocessing...")
pre_df = scLENS.preprocess(ndf)

cur_dev = has_cuda() ? "gpu" : "cpu"
println("Using device: $cur_dev")

println("Running scLENS embedding...")
sclens_embedding = scLENS.sclens(pre_df, device_=cur_dev)

println("Applying UMAP...")
scLENS.apply_umap!(sclens_embedding)

# Save PCA
println("Saving PCA -> $pca_out")
CSV.write(pca_out, sclens_embedding[:pca_n1])

# Save UMAP
println("Saving UMAP -> $umap_out")
CSV.write(umap_out, DataFrame(sclens_embedding[:umap], :auto))

# Save AnnData
println("Saving h5ad -> $h5ad_out")
scLENS.save_anndata(h5ad_out, sclens_embedding)

# Save UMAP plot
println("Saving plot -> $plot_out")
CairoMakie.activate!(type = "png")
panel = scLENS.plot_embedding(sclens_embedding, pre_df.cell)
save(plot_out, panel)

println("Done!")