## Replicating the environment on another machine

1. Download Miniconda for Python 3 and install it.
2. CLone this repository:

```
git clone https://github.com/bodkan/nea-over-time
```

3. Re-create the Python environment:

```
conda env create -f environment.yml
```

4. Install R packages (I used R 3.4.3).

I wasted hours trying to get R reproducible environment to work using `packrat`.
In the end I gave up - it turned out to be impossible to set up a single
working environment for both Linux and macOS, using both CRAN and Bioconductor
repositories.

```
install.packages(c("devtools", "tidyverse"))

source("https://bioconductor.org/biocLite.R")
biocLite(c("BSgenome.Hsapiens.UCSC.hg19", "VariantAnnotation", "rtracklayer",
"GenomicRanges", "IRanges"))

devtools::install_github("bodkan/slimr")
```

## Creating Python environment
```
conda create --name nea-over-time python=3.6
source activate rubber-ducks
```

### Setting up Python environment

```
conda install numpy pandas jupyter jupyterlab
conda install -c bioconda snakemake
conda install -c conda-forge msprime altair vega_datasets
```

Then export the dependency file `environment.yml` using (appnope is a weird
macOS-specific package that does not exist on Linux and is not needed there):

```
conda env export --no-builds | grep -v appnope > environment.yml
```
```
conda install numpy pandas jupyter jupyterlab
conda install -c conda-forge pysam pybedtools
```

```
install.packages("devtools", "tidyverse", "magrittr"))

source("https://bioconductor.org/biocLite.R")
biocLite(c("BSgenome.Hsapiens.UCSC.hg19", "VariantAnnotation", "rtracklayer",
"GenomicRanges", "IRanges"))

devtools::install_github("bodkan/admixr")
devtools::install_github("bodkan/slimr")
```
