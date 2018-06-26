## Replicating the environment on another machine

4. Install R packages (I used R 3.4.3).

I wasted hours trying to get R reproducible environment to work using `packrat`.
In the end I gave up - it turned out to be impossible to set up a single
working environment for both Linux and macOS, using both CRAN and Bioconductor
repositories.

You'll have to install all R dependencies manually like this:

```
install.packages(c("devtools", "tidyverse"))

source("https://bioconductor.org/biocLite.R")
biocLite(c("BSgenome.Hsapiens.UCSC.hg19", "VariantAnnotation", "rtracklayer",
"GenomicRanges", "IRanges"))

devtools::install_github("bodkan/slimr")
devtools::install_github("bodkan/admixr", ref = "v0.1")

install.packages(c('repr', 'IRdisplay', 'evaluate', 'crayon', 'pbdZMQ', 'devtools', 'uuid', 'digest'))
devtools::install_github('IRkernel/IRkernel')
IRkernel::installspec()
```
