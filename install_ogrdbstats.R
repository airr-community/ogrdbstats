# Install ogrdbstats - intended for use in building the Docker image

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
options(repos = BiocManager::repositories())
remotes::install_github("airr-community/ogrdbstats")
