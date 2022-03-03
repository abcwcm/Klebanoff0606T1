[![DOI](https://zenodo.org/badge/465691691.svg)](https://zenodo.org/badge/latestdoi/465691691)

# SIFT-seq data for patient 0606T1

R data package related to 

>Chandran et al. *Immunogenicity and therapeutic targeting of a public neoantigen derived from mutated PIK3CA* Nat Medicine (2022).

This package contains access to the processed TCR-seq and scRNA-seq data of T cells obtained from patient 0606T1.

## How to use it

```
## install
devtools::install_github("abcwcm/Klebanoff0606T1")
```

Upon installation, the processed data, e.g. in the form of `SingleCellExperiment` objects, can be loaded thus:

```
## load SingleCellExperiment objects
sce.0606 <- Klebanoff0606T1::load_0606T1shared()

## load results of differential gene expression comparisons
## see 
Klebanoff0606T1::load_DE_results() # loads an object named `delist.both`
de.0606 <- delist.both; rm(delist.both)
```

For more details, see the [code repository](https://github.com/abcwcm/Chandran2021) and the `Rmd` files detailing how the data underlying the figures were obtained.
