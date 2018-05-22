## Replicating the environment on another machine

1. Download Miniconda for Python 3 and install it.
2. Clone this repository:

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
devtools::install_github("bodkan/admixr", ref = "v0.1")
```

## Creating Python environment
```
conda create --name nea-over-time python=3.6
source activate nea-over-time
```

### Setting up Python environment

```
conda install pandas jupyter jupyterlab matplotlib
conda install -c bioconda pysam pybedtools
```

Then export the dependency file `environment.yml` using (appnope is a weird
macOS-specific package that does not exist on Linux and is not needed there):

```
conda env export --no-builds | grep -v appnope > environment.yml
```
