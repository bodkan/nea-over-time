To reproduce absolutely everything from scratch, You'll need to install
the following dependencies:

## Python

I used Python version 3.6.5 and the following Python modules:

```
pip install numpy pandas msprime pybedtools jupyter
```

## R

I used R version 3.4.3.

Packages from CRAN:
```
install.packages(c("devtools", "broom", "dunn.test", "forcats", "future",
                   "ggbeeswarm", "ggrepel", "here", "latex2exp", "magrittr",
                   "Metrics", "modelr", "parallel", "purrr", "stringr",
                   "tidyverse"))
```

Packages from Bioconductor:
```
source("https://bioconductor.org/biocLite.R")
biocLite(c("biomaRt", "VariantAnnotation", "BSgenome.Hsapiens.UCSC.hg19",
           "GenomicRanges",  "rtracklayer"))
```

Packages from GitHub:
```
devtools::install_github("thomasp85/patchwork")
devtools::install_github("bodkan/bdkn")
devtools::install_github("bodkan/slimr", ref = "v0.1")
devtools::install_github("bodkan/admixr", ref = "v0.1")
```

To be able to run Jupyter notebooks that contain all my analses and figures,
you will also need to install [IRkernel](https://irkernel.github.io).

## SLiM

I used SLiM v2.6. Be aware that SLiM introduced some backwards incompatible
changes since it's 2.0 release, so make sure to use exactly version 2.6.
