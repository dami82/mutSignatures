# mutSignatures
Thank you for your interest in `mutSignatures`! 

- version: latest (dev) version - 2.1.3
- previous versions are available as alternative branches

### Notes about the R library

- integrated R-based computational framework aimed at deciphering DNA mutational signatures

- provides advanced functions for importing DNA variants, computing mutation types, and 
extracting mutational signatures via non-negative matrix factorization (NMF) 

- I/O: accepts multiple types of input data (VCF, MAF), is compatible with non-human genomes, 
and supports the analysis of non-standard mutation types, such as tetra-nucleotide mutation types. 


### Links

- **Vignette on CRAN**: <https://cran.r-project.org/web/packages/mutSignatures/vignettes/get_sarted_with_mutSignatures.html>

- **Peer-reviewed paper (Sci Rep. 2020 Oct 26;10(1):18217)**: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7589488/>

### Installation from GitHub

In R, you can run:
```
devtools::install_github("dami82/mutSignatures", force = TRUE, build_vignettes = TRUE)
```

### Authors
Damiano Fantini, Vania Vidimar, Yanni Yu, Salvatore Condello, Joshua J Meeks
The package was developed at Northwestern University, Chicago, IL (2017-2020)
