# Virusparies

VirusHunterGatherer is a data-driven tool designed for high-throughput virus discovery.
The process involves Virushunter conducting sensitive homology-based detection of
viral sequence reads in unprocessed data. This approach identifies the most conserved regions
of a virus, which serve as seeds in the Virusgatherer step for a seed-based assembly
of full-length viral genome sequences.

The Virusparies package provides a set of plotting functions tailored for visualizing
**VirushunterGatherer Hittables**. The name draws inspiration from the hunter-gatherer metaphor,
with "paries" derived from Latin meaning "wall". It symbolizes the parietal art left by
ancient hunters and gatherers on walls, summarizing their stories and beliefs.

VirusHunterGatherer is available on: https://github.com/lauberlab/VirusHunterGatherer .

- [Installation](https://github.com/SergejRuff/Virusparies/tree/main?tab=readme-ov-file#installation)
- [Overview](https://github.com/SergejRuff/Virusparies/tree/main#overview)
- [Details](https://github.com/SergejRuff/Virusparies/tree/main?tab=readme-ov-file#details)
- [Contributions](https://github.com/SergejRuff/Virusparies/tree/main?tab=readme-ov-file#contributions)


## Installation

You can choose to install the development version of Virusparies from GitHub:

``` r
### First install the "remotes" package

install.packages("remotes")

### Then install the Virusparies package

remotes::install_github("SergejRuff/Virusparies")

```
***!*** Bioconductor version coming soon ...


## Overview

Virusparies includes the following functions:

### Import

- `importVirusTable()`  Import VirusHunterGatherer Hittables into R.


### VirushunterGatherer Plots:

- `VhgBoxplot()` Boxplot plotting RefSeq identity, evalues or contig length for each group.
- `VhgIdenFacetedScatterPlot()` Faceted Scatter Plot for Reference Identity vs -log10 of Reference E-value.
- `VhgIdentityScatterPlot()` Scatter Plot for Reference Identity vs -log10 of Reference E-value.


### Virushunter only Plots:

- `vhRunsBarplot()` Bar plot showing how many unique Runs map against each virus.
- `vhSumHitsBarplot()` Bar plot for the sum of hits for each virus found in group.


### Graphical Tables(gt):

- `vhRunsTable()` Generate a gt-table for vhRunsBarplot.
- `vhgTabularRasa()` Generate custome gt-tables.


### Export:

- `exportVirusGt()`  Export Graphical Tables.
- `exportVirusPlot()` Export plots.


## Details
 
The Virusparies package provides a set of plotting functions tailored for visualizing
VirushunterGatherer Hittables and generating graphical tables (gt) for the results generated by Virusparies or the user.
It also contains functions for import of VirusHunterGatherer Hittables and export of both plots and graphical tables.
Here the functionality of each function shall be explained in detail.

First we load the package into R:

``` r
### load Virusparies into R

library(Virusparies)
```

### Import

VirusHunterGatherer Hittables are tab-separated values (TSV) files, which can be imported with the  `importVirusTable()` function.
Here we will use the example Virushunter and VirusGatherer Hittable files included in  the Virusparies package as an example:

``` r
### Import VirusHunter Hittable.

path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
vh_file <- importVirusTable(path)

print(head(vh_file))  # print head of hunter files

### Import VirusGatherer Hittable.

path2 <- system.file("extdata", "virusgatherer.tsv", package = "Virusparies")
vg_file <- importVirusTable(path2)

print(head(vg_file))  # print head of gatherer files
```

### VirusHunterGatherer Plot - VhgBoxplot

``` r
### Load VirusHunter File

path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
vh_file <- importVirusTable(path)

### Plot 1 for evalues
plot1 <- VhgBoxplot(vh_file, x_column = "best_query", y_column = "ViralRefSeq_E")
plot1

### Plot 2 for identity
plot2 <- VhgBoxplot(vh_file, x_column = "best_query", y_column = "ViralRefSeq_ident")
plot2

### Plot 3 custom arguments used
plot3 <- VhgBoxplot(vh_file,
                   x_column = "best_query",
                   y_column = "ViralRefSeq_E",
                   theme_choice = "grey",
                   subtitle = "Custom subtitle: Identity for custom query",
                   xlabel = "Custom x-axis label: Custom query",
                   ylabel = "Custom y-axis label: Viral Reference Evalue in -log10 scale",
                   legend_position = "right")
plot3

### Import gatherer files
path2 <- system.file("extdata", "virusgatherer.tsv", package = "Virusparies")
vg_file <- importVirusTable(path2)

### Plot 4: virusgatherer plot with SRA_run as custom grouping
plot4 <- VhgBoxplot(vg_file,x_column = "SRA_run",y_column = "ViralRefSeq_E")
plot4

### Plot 5: Virusgatherer plot for SRA_runs agains contig length
plot5 <- VhgBoxplot(vg_file,x_column = "SRA_run",y_column = "contig_len")
plot5
```

### VirusHunterGatherer Plot - VhgIdenFacetedScatterPlot

``` r

### Load VirusHunter File

path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
vh_file <- importVirusTable(path)

### Generate Plot

plot <- VhgIdenFacetedScatterPlot(vh_file,cutoff = 1e-5)
plot

```


### VirusHunterGatherer Plot - VhgIdentityScatterPlot

``` r

### Load VirusHunter File

path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
vh_file <- importVirusTable(path)

### Basic plot

plot <- VhgIdentityScatterPlot(vh_file,cutoff = 1e-5)
plot(plot)

```


### VirusHunter Plot - vhRunsBarplot

``` r

### Load VirusHunter File

path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
vh_file <- importVirusTable(path)

### Generate Plot

plot <- vhRunsBarplot(vh_file,cut = 1e-5)
plot


```


### VirusHunter Plot - vhSumHitsBarplot

``` r

### Load VirusHunter File

path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
vh_file <- importVirusTable(path)

### Generate Plot

plot <- vhSumHitsBarplot(vh_file,cut = 1e-5)
plot

```


### VirusHunter GT - vhRunsTable

`vhRunsTable()` takes VirusHunter Files as Input and generates a graphical table providing the user with information about which Run has found which virus group.
This makes it a nice complement to the `vhRunsBarplot()`function as they provide the user with the same information but with more detail about the exact number of unique Runs in the plot
and the Run-ID in the table.

``` r

### Load VirusHunter File

path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
vh_file <- importVirusTable(path)

### Generate table with defaul arguments

table <- vhRunsTable(vh_file,cut = 1e-5)
table

```


### GT - vhgTabularRasa

This function creates a formatted table using the `gt` package, based on input data with specified column names.
It is particularly useful for generating tables that cannot be produced with `vhRunsTable`.

``` r

### Load VirusHunter File

path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
vh_file <- importVirusTable(path)

### Plot boxplot for "identity"

identity <- VhgBoxplot(vh_file,y_column = "ViralRefSeq_ident")

# Generate table

vhgTabularRasa(identity$summary_stats)


```

### Export

#### Export Plots

Plots can be exported in different formats via `exportVirusPlot()`.
Supported devices include "eps", "ps", "tex", "pdf", "jpeg", "tiff", "png", "bmp", "svg", or "wmf" (Windows only).
When 'device' argument is set to NULL, the file extension in 'filename' is used to determine the device.

``` r
### Load VirusHunter File

path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
vh_file <- importVirusTable(path)

### Generate Basic plot

plot <- VhgIdentityScatterPlot(vh_file,cutoff = 1e-5)


### Export plot

exportVirusPlot(plot=plot,file_name="testplot.png",width=8,height=6,units="in")


```
***!*** Depending on the plot, the final image might be cropped or truncated.
We recommend experimenting with height, width, and resolution.

#### Export Graphical Tables

`exportVirusPlot()` implements the `gtsave` function from the gt-package to export graphical tables in different formats.
This feature is in an experimental phase and may not currently function as expected. 
***!*** Google Chrome or a Chromium browser is required for png and pdf file export.

``` r

### Load VirusHunter File

path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
vh_file <- importVirusTable(path)

### Using first 10 rows of SRA_run,num_hits,bestquery,ViralRefSeq_E and Identity col.

vh_file_part <- vh_file[c(1:10),c(1,7,9,10,11)]

### Generating a gt

table <- vhgTabularRasa(vh_file_part,title = "first 10 rows of vh_file",subtit =
"example for any table",names_ = c("Runs","Number of Contigs","Best Query Result",
"Reference E-Value","Refrence Identity"))

### Export gt as docx file

exportVirusGt(gtable=table,filename="vh_parttable.docx")

```

## Contributions

Sergej Ruff formulated the idea behind Virusparies and was responsible for its implementation.

Chris Lauber and Li Chuin Chong from Twincore - Centre for Experimental and Clinical Infection Research provided ideas for improvements. 
Chris Lauber is the main developer behind the VirusHuntergatherer software and the Group Leader of the Computational Virology working group at the Institute for Experimental Virology, TWINCORE.
