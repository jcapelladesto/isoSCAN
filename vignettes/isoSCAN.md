---
title: "isoSCAN Package Vignette"
author: "Jordi Capellades"
date: "2020-07-01"
output:
  html_vignette:
    keep_md: yes
vignette: >
  %\VignetteIndexEntry{isoSCAN Package Vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction
This vignette contains a basic walkthrough of the functionalities of the `isoSCAN` package. The package is designed to automatically extract the abundances of isotopologues of a targeted list of compounds. It is capable of doing so in both low- and high-resolution data, though depending on the resolution the requirements for the input are different.

### Targeted compound list format
This package requires a specific __targeted compound list__ format that will be used in _autoQ_ _formulaTable_ argument. This file can be created on Excel or similar software and then imported into R via `read.csv`. _formulaTable_ __must__ contain the following column names in no specific order:

* __CompoundName__ is the name of the compound or metabolite quantified

* __mz__ of the monoisotopic ion

* __RT__ retention time value in seconds

* __Formula__ of the compound. __NOTE:__ This formula must match the derivatized _Formula_ including derivatization modifications in the case of high-resolution data.

* __NumAtoms__ determining the number of compounds to be quantified

* Other columns will be ignored 


### Other isoSCAN functionalities
The rest of the functions are used for processing the raw data for either quantification or plotting. The package currently contains the following functions.

* Quantifcation:

    * `autoQ`

* Data transformation

    * `QTransform`
    
    * `simplifybyPpm`
        
    * `sumIsotopologues`
    
* Plotting

    * `rawPlot`

    * `meanRawPlot`

    * `metBarPlot`

### Creating the formulaTable
Before starting with file processing, we need to load the _targeted compounds_ as a _formulaTable_ data frame. This can be done either with `read.table` or `read.csv` functions. Make sure that the file contains it contains the columns as listed in the section above.

The package includes examples for both Low- and High-Resolution:


```r
library(isoSCAN)
data("formulaTables")

# Low-Resolution (e.g. nominal mass accuracy)
formulaTable_lowres <- formulaTables[which(formulaTables$Instrument=="Quadrupole"),]

formulaTable_lowres
##   CompoundName    RT    mz  Formula NumAtoms Instrument
## 1         ILeu 357.0 276.0 C6H13NO2        6 Quadrupole
## 3          Leu 373.2 276.0 C6H13NO2        6 Quadrupole
## 4          Gly 382.2 292.0  C2H5NO2        2 Quadrupole
## 5    Succinate 384.0 263.1   C4H4O4        4 Quadrupole
## 6     Fumarate 406.8 261.1   C4H2O4        4 Quadrupole
## 7          Ala 418.2 306.0  C3H7NO2        3 Quadrupole
## 8          Ser 421.8 322.0  C3H7NO3        3 Quadrupole
```



```r
# High-Resolution (Orbitrap, or qTOF)
formulaTable_orbi <- formulaTables[which(formulaTables$Instrument=="Orbitrap"),]

formulaTable_orbi
##    CompoundName    RT       mz      Formula NumAtoms Instrument
## 9           Val 354.6 262.1653 C11H28NO2Si2        5   Orbitrap
## 11         ILeu 396.6 276.1809 C12H30NO2Si2        6   Orbitrap
## 12          Leu 411.6 276.1809 C12H30NO2Si2        6   Orbitrap
## 13          Gly 421.8 292.1579 C11H30NO2Si3        2   Orbitrap
## 14    Succinate 426.6 263.1129  C10H23O4Si2        4   Orbitrap
## 15     Fumarate 453.0 261.0973  C10H21O4Si2        4   Orbitrap
## 16          Ser 459.0 322.1684 C12H32NO3Si3        3   Orbitrap
## 17          Ala 460.2 306.1735 C12H32NO2Si3        3   Orbitrap
```


### Creating and loading mz(X)ML files
The first step is file format transformation, `isoSCAN` uses `mzR` package in order to read MS files. Therefore, you will have to transform the raw data from vendor format into __mz(X)ML__ format using __Proteowizard MSconvert__ (or similar tools), so they can be read by the `mzR` R package.
There is an important parameter to consider in MSconvert depending on the nature of the data resolution:
* In the case of Low-resolution. Transform the data mantaining __profile format__. This is essential for peak quantification. (e.g. _peakPicking=False_ in MSconvert)
* In the case of High-resolution, __please use centroiding__ (e.g. _peakPicking= True_ in MSconvert)

Then, we need to locate the folder in which these files are found and list them in a vector.


```r
#run this: setwd("./mydatafolder")
SampleFiles  <- list.files(pattern="\\.mz(X)?ML")
```

This package also includes sample mzML data files to be used for testing:

```r
# Low-resolution files
SampleFiles_lowres <- list.files(system.file("extdata",package = "isoSCAN"),
                                 full.names = T,pattern = "lowres")

#High-resolution files
SampleFiles_orbi <- list.files(system.file("extdata",package = "isoSCAN"),
                               full.names = T,pattern = "orbi")
```

### Processing files
Now we can call `autoQ` function that will process the files and look for the isotopologues for each compound found in the `formulaTable`.
Additionally, other parameters need to be indicated as stated in _help(autoQ)_. This parameters refer to peak width and number of scans recorded, together with signal-to-noise ratio and mass error.

#### Low-resolution data
In the case of low-resolution data. Please remember to use them in _Profile_ format as it eases the process of peak finding. 



```r
head(integrations)
##   CompoundName      m.z Isotopologue maxo_lowres1_C12 area_lowres1_C12
## 1         ILeu 276.0000          M+0       27803.3535       235102.860
## 2         ILeu 277.0034          M+1        7990.1074       107147.079
## 3         ILeu 278.0067          M+2               NA               NA
## 4         ILeu 279.0101          M+3        1756.6731               NA
## 5         ILeu 280.0134          M+4         968.5995         9400.776
## 6         ILeu 281.0168          M+5         251.8362               NA
##   maxo_lowres1_C13 area_lowres1_C13 maxo_lowres2_C12 area_lowres2_C12
## 1       63631.6445        566743.25       44112.3359        373584.46
## 2       17535.6797        277193.54       12711.3623        185009.87
## 3        9129.6455        126876.92        7473.7661         66831.24
## 4        2112.3086         24435.29               NA               NA
## 5        1026.8778               NA         889.4612               NA
## 6         407.8593               NA         392.3322               NA
##   maxo_lowres2_C13 area_lowres2_C13
## 1       63551.7461         546198.4
## 2       16982.4316         276727.2
## 3        6969.1235         122360.4
## 4        1641.4408               NA
## 5        1131.9938               NA
## 6         450.8288               NA
```

#### High-resolution data
###### TODO: explain enviPat

```r
data(isotopes, package="enviPat")
isotopes[isotopes$isotope=="13C",] # both rows required
##     element isotope     mass abundance ratioC
## 11        C     13C 13.00335    0.0107      0
## 297   [13]C     13C 13.00336    1.0000      0
```





```r
head(integrations)
##   CompoundName         m.z      abundance Isotopologue          ppm
## 1          Val 262.1653085 100.0000000000          M+0 0.4853617595
## 2          Val 263.1647635  10.0055163842          M+0 3.1890036771
## 3          Val 263.1688703  12.8262035042          M+1           NA
## 4          Val 264.1622733   7.0433925761          M+0 1.0826693676
## 5          Val 264.1687360   1.6578098230          M+2 1.4853580079
## 6          Val 265.1657067   0.8633555587          M+1 1.8556071808
##   maxo_orbi1_C12 area_orbi1_C12 ppm maxo_orbi1_C13 area_orbi1_C13          ppm
## 1  8976119.00000  60551262.4098  NA             NA             NA 0.4458851181
## 2  1770118.37500             NA  NA             NA             NA           NA
## 3             NA             NA  NA             NA             NA 1.9182257184
## 4   594159.75000   3854748.8445  NA             NA             NA 0.8516176149
## 5    95426.71875    575899.2997  NA             NA             NA 1.8319271582
## 6    59326.93750    334466.2617  NA             NA             NA 1.2801636496
##   maxo_orbi2_C12 area_orbi2_C12          ppm maxo_orbi2_C13 area_orbi2_C13
## 1 1419617.625000             NA 0.3689558998  9450046.00000  65156489.5397
## 2             NA             NA 2.6091848267  1751672.25000             NA
## 3  113213.656250             NA 1.6863017975   126260.34375             NA
## 4   91571.453125             NA 0.7938546767   625612.06250   4162402.1243
## 5   14997.597656             NA 1.4853580079    94064.78906    609487.0947
## 6    8817.168945             NA 2.3159620059    62637.26562    366316.0257
```

This processes each file independently, looking for "good-shape" peaks and obtaining both the area and max intensity scan (Maxo) for each isotopologue, if the area cannot be calculated (due to noise or peak shape) then only the Maxo is returned.

Once finished, we can plot them or transform the values for exportation. 

## Plotting results and value transformation
The `metBarPlot` function is designed to plot values in a barplot including standard deviation error bars. The arguments for value are the following:

* _groups_ Sample groups. This should be a vector with the groups of experiments, matching the same order in `list.files(pattern=".mzXML")`. This vector can be created using `gsub` function and others. See the example:

* _val.to.plot_ Values to use for plotting. Either "Area" or "Maxo".

* _ylabel_ The text that should be shown in the Y axis (i.e. Intensity, Area). Otherwise, _val.to.plot_ value will be used. 


```r
mygroups <- isoSCAN:::rmfileExt(SampleFiles_orbi,"\\.mzML")
mygroups <- gsub(".*_","",mygroups)
```



```r
metBarPlot(autoQres=integrations, groups = mygroups, val.to.plot="area")
```


```r
# Example of a single compound
integrations <- sumIsotopologues(integrations) #only for high-res
metBarPlot(autoQres=integrations[which(integrations$CompoundName=="Gly"),],
           groups = mygroups, val.to.plot="area")
## Warning in metBarPlot(autoQres = integrations[which(integrations$CompoundName
## == : 'ylabel' value is NULL. Using 'val.to.plot' value as default label for y
## axis.
## Not quantified
## Not quantified
## Not quantified
## Not quantified
```


```r
#if you want to save the plots to a file:
metBarPlot(autoQres=integrations, groups = mygroups, val.to.plot="area",
					 topdf="./metBarPlot_results.pdf",height=10,width=18,pointsize=16)
#modify height, width and pointsize accordingly to fit your output size
```

`metBarPlot` can also digest the data frame produced by `QTransform`. Read the function help for further information. 


```r
trans_integ <- QTransform(integrations,val.trans = "P")
metBarPlot(autoQres=trans_integ, groups = mygroups, val.to.plot="area")
```
##### TODO: ERROR?

## Extra: Raw data plotting
`rawPlot` and `meanRawPlot` functions should be used for quality control purposes. They are useful to check for moving peaks, noisy spots or saturated peaks.
The first function will print a "spectra heatmap", m/z and retention time in the x,y-axis respectively, with points coloured depending on the scan intensity value. 


####### TODO : SOLVE  could not find function "rmfileExt"

While the `rawPlot` function will run through each file and plot the raw spectra for all compounds in a single pdf file, `meanRawPlot` will calculate the average spectra for all files and generate a single plot for each compound.


```r
# Example plot for the first compound, note that now we use ALL files
meanRawPlot(SampleFiles=SampleFiles_lowres[c(1,3)],
            formulaTable=formulaTable_lowres[1:2,], 
            RTwin=5)

```
####### TODO : SOLVE  could not find h

#### Saving plots into PDF documents
Both functions contain the `topdf` argument, that is encouraged to be used in order to save the files into _PDF format_ instead of being shown in the R plotting default device. In the case of `meanRawPlot` indicate the name of the output PDF file desired (`topdf=/plot_folder/mean_raw_spectra.pdf`), whereas in the case of `rawPlot` should be used as: `topdf="C:/Users/User/Desktop/plot_folder` and all plots generated will be saved into *plot_folder* mantaining the sample name for the pdf file name.
