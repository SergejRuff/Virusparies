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
    - [Import](https://github.com/SergejRuff/Virusparies/blob/main/README.md#import-1)
    - [VirusHunterGatherer Plot - VhgBoxplot](https://github.com/SergejRuff/Virusparies/blob/main/README.md#virushuntergatherer-plot---vhgboxplot)
    - [VirusHunterGatherer Plot - VhgIdenFacetedScatterPlot](https://github.com/SergejRuff/Virusparies/blob/main/README.md#virushuntergatherer-plot---vhgidenfacetedscatterplot)
    - [VirusHunterGatherer Plot - VhgIdentityScatterPlot](https://github.com/SergejRuff/Virusparies/blob/main/README.md#virushuntergatherer-plot---vhgidentityscatterplot)
    - [VirusHunter Plot - vhRunsBarplot](https://github.com/SergejRuff/Virusparies/blob/main/README.md#virushunter-plot---vhrunsbarplot)
    - [VirusHunter Plot - vhSumHitsBarplot](https://github.com/SergejRuff/Virusparies/blob/main/README.md#virushunter-plot---vhsumhitsbarplot)
- [Citation](https://github.com/SergejRuff/Virusparies/blob/main/README.md#citation)
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
 
The Virusparies package provides a set of plotting functions tailored for visualizing VirushunterGatherer hittables and generating graphical tables (gt) for the results generated by Virusparies or the user. It also contains functions for importing VirusHunterGatherer hittables and exporting both plots and graphical tables. Here, the functionality of each function will be explained in detail.

First, we load the package into R:

``` r
### load Virusparies into R

library(Virusparies)
```

### Import

VirusHunterGatherer hittables are tab-separated values (TSV) files, which can be imported with the `importVirusTable()` function. Here, we will use the example VirusHunter and VirusGatherer hittable files included in the Virusparies package as an example:

``` r
### Import VirusHunter Hittable.

path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
vh_file <- importVirusTable(path)

print(head(vh_file))  # print head of hunter files

#>      SRA_run SRA_sample SRA_study   host_taxon host_taxid num_queries num_hits  best_E
#> 1 SRR10822543 SRS5937532 SRP239389 Homo_sapiens       9606           8       12 4.5e-05
#> 2 SRR12567985 SRS7305728 SRP279687 Homo_sapiens       9606           6        9 9.7e-13
#> 3 SRR12567985 SRS7305728 SRP279687 Homo_sapiens       9606           6        5 1.4e-08
#> 4 SRR10822560 SRS5937549 SRP239389 Homo_sapiens       9606          11        3 1.1e-08
#> 5 SRR12567985 SRS7305728 SRP279687 Homo_sapiens       9606           6        8 3.2e-10
#> 6 SRR10822546 SRS5937535 SRP239389 Homo_sapiens       9606           7       20 1.0e-01
#>       best_query ViralRefSeq_E ViralRefSeq_ident ViralRefSeq_aLen.sLen ViralRefSeq_contigs
#> 1 Anello_ORF1core      2.59e-38              93.8             192 / 194                   1
#> 2 Anello_ORF1core      8.44e-36              77.3             225 / 227                   6
#> 3 Anello_ORF1core      3.95e-34              70.1             231 / 232                   1
#> 4 Anello_ORF1core      3.69e-31              68.6             210 / 216                   1
#> 5 Anello_ORF1core      5.85e-31              76.9             195 / 201                   3
#> 6 Anello_ORF1core      8.63e-31              79.4             189 / 237                   1
#>                                                                               ViralRefSeq_subject
#> 1 gi:1464307144|Torque teno mini virus 11 isolate LIL-y2 ORF2, ORF1, and ORF3 genes, complete cds
#> 2                   gi:1464307216|Torque teno midi virus 11 DNA, complete genome, isolate: MDJN47
#> 3 gi:1464307144|Torque teno mini virus 11 isolate LIL-y2 ORF2, ORF1, and ORF3 genes, complete cds
#> 4                                          gi:134133206|Torque teno midi virus 1, complete genome
#> 5                   gi:1464307186|Torque teno midi virus 5 DNA, complete genome, isolate: MDJHem2
#> 6                   gi:1464307221|Torque teno midi virus 12 DNA, complete genome, isolate: MDJN51
#>                           ViralRefSeq_taxonomy date_analyzed
#> 1  taxid:2065037|Betatorquevirus|Anelloviridae    2022-05-17
#> 2 taxid:2065052|Gammatorquevirus|Anelloviridae    2022-05-06
#> 3  taxid:2065037|Betatorquevirus|Anelloviridae    2022-05-06
#> 4  taxid:687379|Gammatorquevirus|Anelloviridae    2022-05-17
#> 5 taxid:2065046|Gammatorquevirus|Anelloviridae    2022-05-06
#> 6 taxid:2065053|Gammatorquevirus|Anelloviridae    2022-05-17


### Import VirusGatherer Hittable.

path2 <- system.file("extdata", "virusgatherer.tsv", package = "Virusparies")
vg_file <- importVirusTable(path2)

print(head(vg_file))  # print head of gatherer files

#>       SRA_run SRA_sample SRA_study   host_taxon host_taxid                 contig_id
#>  1   ERR206007  ERS074208 ERP000373 Homo sapiens       9606   ERR206007_cap3_Contig-1
#>  2   ERR206007  ERS074208 ERP000373 Homo sapiens       9606   ERR206007_cap3_Contig-2
#>  3   ERR206007  ERS074208 ERP000373 Homo sapiens       9606   ERR206007_cap3_Contig-3
#>  4   ERR206007  ERS074208 ERP000373 Homo sapiens       9606   ERR206007_cap3_Contig-4
#>  5   ERR206021  ERS074222 ERP000373 Homo sapiens       9606   ERR206021_cap3_Contig-1
#>  6 SRR10822543 SRS5937532 SRP239389 Homo sapiens       9606 SRR10822543_cap3_Contig-1
#>    contig_len ViralRefSeq_E ViralRefSeq_ident ViralRefSeq_aLen
#>  1        603      4.22e-65            92.453              106
#>  2        461      3.22e-65            85.612              139
#>  3        364      2.46e-53            76.531               98
#>  4        334      9.87e-67            93.636              110
#>  5        323      8.91e-45            65.421              107
#>  6       3321      0.00e+00            86.058              789
#>                                                       ViralRefSeq_subject
#>  1                            acc:YP_009505712|Orf1 [Torque teno virus 5]
#>  2 acc:YP_003587853|hypothetical protein TTV10_gp4 [Torque teno virus 10]
#>  3                            acc:YP_003587868|ORF1 [Torque teno virus 3]
#>  4                           acc:YP_009505715|Orf1 [Torque teno virus 11]
#>  5        acc:YP_009505729|unnamed protein product [Torque teno virus 24]
#>  6                        acc:YP_009173866|polymerase [Hepatitis B virus]
#>                                                                                                ViralRefSeq_taxonomy
#>  1                                              taxid:687344|Alphatorquevirus homin5|Alphatorquevirus|Anelloviridae
#>  2                                             taxid:687349|Alphatorquevirus homin10|Alphatorquevirus|Anelloviridae
#>  3                                              taxid:687342|Alphatorquevirus homin3|Alphatorquevirus|Anelloviridae
#>  4                                                                      taxid:687350|Alphatorquevirus|Anelloviridae
#>  5                                             taxid:687363|Alphatorquevirus homin24|Alphatorquevirus|Anelloviridae
#>  6 taxid:10407|Orthohepadnavirus|Hepadnaviridae|Blubervirales|Revtraviricetes|Artverviricota|Pararnavirae|Riboviria
#>    date_analyzed
#>  1    2024-05-18
#>  2    2024-05-18
#>  3    2024-05-18
#>  4    2024-05-18
#>  5    2024-05-17
#>  6    2024-05-18


```

### VirusHunterGatherer Plot - VhgBoxplot

#### Boxplot 1: "ViralRefSeq_E"

``` r


### Plot 1 for evalues
plot1 <- VhgBoxplot(vh_file, x_column = "best_query", y_column = "ViralRefSeq_E")
plot1


```
#### Boxplot 2: "ViralRefSeq_ident"

``` r

### Plot 2 for identity
plot2 <- VhgBoxplot(vh_file, x_column = "best_query", y_column = "ViralRefSeq_ident")
plot2

``` 
#### Boxplot 3: Customization

``` r

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

```

#### Boxplot 4: "SRA_run" as custom group

``` r
### Plot 4: virusgatherer plot with SRA_run as custom grouping
plot4 <- VhgBoxplot(vg_file,x_column = "SRA_run",y_column = "ViralRefSeq_E")
plot4

```

#### Boxplot 5: "contig_len" (Gatherer Tables only)

``` r

### Plot 5: Virusgatherer plot for SRA_runs agains contig length
plot5 <- VhgBoxplot(vg_file,x_column = "SRA_run",y_column = "contig_len")
plot5

``` 


### VirusHunterGatherer Plot - VhgIdenFacetedScatterPlot

``` r

### Generate Plot

plot <- VhgIdenFacetedScatterPlot(vh_file,cutoff = 1e-5)
plot

```


### VirusHunterGatherer Plot - VhgIdentityScatterPlot

``` r

### Basic plot

plot <- VhgIdentityScatterPlot(vh_file,cutoff = 1e-5)
plot(plot)

```


### VirusHunter Plot - vhRunsBarplot

``` r

### Generate Plot

plot <- vhRunsBarplot(vh_file,cut = 1e-5)
plot


```


### VirusHunter Plot - vhSumHitsBarplot

``` r

### Generate Plot

plot <- vhSumHitsBarplot(vh_file,cut = 1e-5)
plot

```


### VirusHunter GT - vhRunsTable

The `vhRunsTable()` function takes VirusHunter files as input and generates a graphical table, providing the user with information about which run has found which virus group. This makes it a valuable complement to the `vhRunsBarplot()` function. While `vhRunsBarplot()` provides the information in a plot that quantifies the number of unique runs finding a virus group, `vhRunsTable()` presents the same information in table form, showing which runs are found along with their names (SRA accessions, FASTQ).


``` r


### Generate table with defaul arguments

table <- vhRunsTable(vh_file,cut = 1e-5)
table

```


### GT - vhgTabularRasa


This function creates a formatted table using the gt package, based on input data with specified column names. It is particularly useful for generating tables that cannot be produced with `vhRunsTable`.


``` r


### Plot boxplot for "identity"

identity <- VhgBoxplot(vh_file,y_column = "ViralRefSeq_ident")

# Generate table

vhgTabularRasa(identity$summary_stats)


```

### Export

#### Export Plots

Plots can be exported in various formats via the `exportVirusPlot()` function. Supported formats include "eps", "ps", "tex", "pdf", "jpeg", "tiff", "png", "bmp", "svg", and "wmf" (Windows only). When the device argument is set to NULL, the file extension in filename is used to determine the export format.


``` r


### Generate Basic plot

plot <- VhgIdentityScatterPlot(vh_file,cutoff = 1e-5)


### Export plot

exportVirusPlot(plot=plot,file_name="testplot.png",width=8,height=6,units="in")


```
***!*** Depending on the plot, the final image might be cropped or truncated.
We recommend experimenting with height, width, and resolution.

#### Export Graphical Tables

The `exportVirusPlot()` function utilizes the `gtsave` function from the gt package to export graphical tables in various formats. This feature is currently in an experimental phase and may not operate as expected. ***!*** Please note that exporting PNG and PDF files requires Google Chrome or a Chromium-based browser.


``` r


### Using first 10 rows of SRA_run,num_hits,bestquery,ViralRefSeq_E and Identity col.

vh_file_part <- vh_file[c(1:10),c(1,7,9,10,11)]

### Generating a gt

table <- vhgTabularRasa(vh_file_part,title = "first 10 rows of vh_file",subtit =
"example for any table",names_ = c("Runs","Number of Contigs","Best Query Result",
"Reference E-Value","Refrence Identity"))

### Export gt as docx file

exportVirusGt(gtable=table,filename="vh_parttable.docx")

```

## Citation

When utilizing Virusparies in your research or software development, kindly reference the R package using the citation obtained from the `citation()` function:


``` r

### Citation function

citation("Virusparies")

#> To cite package ‘Virusparies’ in publications use:
#> 
#>   Ruff S (2024). _Virusparies: Data Visualisations for
#>   VirusHunterGatherer hittables output_. R package version
#>   1.0.0, <https://github.com/SergejRuff/Virusparies>.

#> A BibTeX entry for LaTeX users is

#>   @Manual{,
#>     title = {Virusparies: Data Visualisations for VirusHunterGatherer hittables output},
#>     author = {Sergej Ruff},
#>     url = {https://github.com/SergejRuff/Virusparies},
#>     year = {2024},
#>     note = {R package version 1.0.0},
#>   }

```

## Contributions

Sergej Ruff formulated the idea behind Virusparies and was responsible for its implementation.

Chris Lauber and Li Chuin Chong from Twincore - Centre for Experimental and Clinical Infection Research provided ideas for improvements. 
Chris Lauber is the main developer behind the VirusHuntergatherer software and the Group Leader of the Computational Virology working group at the Institute for Experimental Virology, TWINCORE.
