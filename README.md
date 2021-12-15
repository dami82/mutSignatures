# mutSignatures

Thank you for your interest in `mutSignatures`! 

- version: latest (dev) version, `2.1.4`
- previous versions are available as alternative branches

### Notes about the R library

- integrated R-based computational framework aimed at deciphering DNA mutational signatures

- provides advanced functions for importing DNA variants, computing mutation types, and 
extracting mutational signatures via non-negative matrix factorization (NMF) 

- I/O: accepts multiple types of input data (VCF, MAF), is compatible with non-human genomes, 
and supports the analysis of non-standard mutation types, such as tetra-nucleotide mutation types. 

### Important differences with previous versions

- the latest version of `mutSignatures` (available here on *GitHub*) fixes a known compatibility issue between
older versions of `mutSignatures` and more recent versions of BSgenome objects. If you get an **error** while 
running the `attachContext()` function, please re-install the latest version of `mutSignatures` (from *GitHub*, see below!) 
and try again.

- compared to mutSignatures version 1.3.1-7, the latest `mutSignature` version replaced the `plot()` method with 
the `msigPlot()` method. Please, make sure to use `mutSignatures::msigPlot()` to build plots (as shown in the
vignette attached to the package) especially when reproducing examples from older vignettes. Thanks!


### Links

- **Vignette on CRAN**: <https://cran.r-project.org/web/packages/mutSignatures/vignettes/get_sarted_with_mutSignatures.html>

- **Peer-reviewed paper** *(Sci Rep. 2020 Oct 26;10(1):18217)*: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7589488/>

### Installation from GitHub

In R, you can run:
```
devtools::install_github("dami82/mutSignatures", force = TRUE, build_vignettes = TRUE)
```

### Authors

- Damiano Fantini, Vania Vidimar, Yanni Yu, Salvatore Condello, Joshua J Meeks
- The package was developed at Northwestern University, Chicago, IL (Department of Urology, 2017-2020)
- The package is still maintained by Damiano Fantini in his free time (Damiano has left academia). Therefore, please
be patient when asking for help of submitting new issues! Thanks.

