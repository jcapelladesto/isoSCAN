---
title: "isoSCAN Package Vignette"
author: "Jordi Capellades"
date: "2017-12-10"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{isoSCAN Package Vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction
This vignette contains a basic walkthrough of the functionalities of the `isoSCAN` package.

This package requires a specific input format for the _compound table_, in which we __must__ find:

* __CompoundName__

* __mz__ of the derivate monoisotopic peak

* __RT__ in seconds

* __NumAtoms__ determining the number of compounds to be quantified

* Other columns will be ignored

The package includes the `buildFormulaTable` function to extract the __NumAtoms__ value from a column named __Formula__.

The rest of the functions are used for processing the raw data for either quantification or plotting. The package currently contains the following functions:

* `autoQ`

* `metBarPlot`

* `QTransform`

* `rawPlot`

* `meanRawPlot`


## Creating the FormulaTable
Before starting with the file processing, we need to load the _compound table_ into the _FormulaTable_ data frame. This can be done either with `read.table` or `read.csv` functions. Make sure that the file contains it contains the columns as listed above.


```r
library(isoSCAN)
formulaTable <- read.table(file=system.file("extdata/formulatable_example.txt", package="isoSCAN"),header=T,sep="\t")

formulaTable$CompoundName <- factor(formulaTable$compound,levels=unique(formulaTable$compound)) # in factor format
formulaTable$RT <- formulaTable$tR*60 # retention time was in minutes
formulaTable$mz <- formulaTable$m.z
```

This table contains a column called `Formula`, we will now extract the number of atoms of the compound using `buildFormulaTable`.


```r
formulaTable$Formula[1:5]
formulaTable <- buildFormulaTable(formulaTable)
formulaTable$NumAtoms[1:5]
```

We now have our formulaTable ready to start the file processing.


```r
formulaTable[1:5,c("CompoundName","NumAtoms", "mz","RT")]
```

## Processing mzXML files
First, the user will have to transform the raw data from vendor format into mz(X)ML format using Proteowizard (or similar tools), so they can be read by the `mzR` R package.
Then, we need to locate the folder in which the mzXML files are found and list them in a vector.


```r
#run this: setwd("./mydatafolder")
lf <- list.files(pattern=".mzXML")
```

Now we can call `autoQ` function that will process the files and look for the isotopologues for each compound found in the `formulaTable`.


```r
integrations <- autoQ(SampleFiles=lf, formulaTable=formulaTable, IsoMass = 1.003355, RTwin = 5, toplot = F, ScanE = 0.1, IntE= 0.4, fit.p = 0.05)
```
This function processes each file independently, looking for "good" peaks and obtaining both the area and max intensity scan (Maxo) for each isotopologue, if the area cannot be calculated, due to noise or peak shape, then only the Maxo is returned.

Once finished, we can plot them or transform the values for exportation. 

## Plotting results and value transformation
The `metBarPlot` function is designed to plot values in a barplot including standard deviation error bars. The arguments for value are the following:

* _groups_ Sample groups. This should be a vector with the groups of experiments, matching the same order in `list.files(pattern=".mzXML")`. This vector can be created using `gsub` function and others. See the example:


```r
mygroups <- lf
mygroups <- gsub("_..?mzXML$","",mygroups) # here we delete all characters after the last "_" until the end of the sample name (check regular expressions in R)
```

* _val.to.plot_ Values to use for plotting. Either "Area" or "Maxo".


```r
metBarPlot(x=integrations, groups = mygroups, val.to.plot="Area")
#if you want to save the plots to a file:
metBarPlot(x=integrations, groups = mygroups, val.to.plot="Area",
					 topdf="./metBarPlot_results.pdf",height=10,width=18,pointsize=16)
#modify height, width and pointsize accordingly to fit your output size

```

* _ylabel_ The text that should be shown in the Y axis (i.e. Intensity, Area). Otherwise, _val.to.plot_ value will be used. 

`metBarPlot` can also digest the data frame produced by `Q2Transform`. Read the function help for further information. 

## Extra: Raw data plotting
`rawPlot` and `meanRawPlot` functions should be used for quality control purposes. They are useful to check for moving peaks, noisy spots or saturated peaks.
The first function will print a "spectra heatmap", m/z and retention time in the x,y-axis respectively, with points coloured depending on the scan intensity value. 


```r
# Example plot from first file and first compound
rawPlot(SampleFiles=lf[1],formulaTable=formulaTable[1,],RTwin=5)
```

While the `rawPlot` function will run through each file and plot the raw spectra for all compounds in a single pdf file, `meanRawPlot` will calculate the average spectra for all files and generate a single plot for each compound.


```r
# Example plot for the first compound, note that now we use ALL files
meanRawPlot(SampleFiles=lf,formulaTable=formulaTable[1,],RTwin=5)
```

Both functions contain the `topdf` argument, that is encouraged to be used in order to save the files into _PDF format_ instead of being shown in the R plotting default device. In the case of `meanRawPlot` indicate the name of the output PDF file desired (`topdf=/plot_folder/mean_raw_spectra.pdf`), whereas in the case of `rawPlot` should be used as: `topdf="C:/Users/User/Desktop/plot_folder` and all plots generated will be saved into *plot_folder* mantaining the sample name for the pdf file name.