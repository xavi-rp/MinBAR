# MinBAR

MinBAR is an R package that aims at (1) defining what is the minimum or optimal background extent necessary to fit good partial SDMs and/or (2) determining if the background area used to fit a partial SDM is reliable enough to extract ecologically relevant conclusions from it.

The idea is to sequentially fit several concentric SDMs, with different diameter each (i.e. buffers), from the centre of the species distribution to the periphery, until a satisfactory model is reached.

The main function of MinBAR is *minba*. In the version 1.0.0 of the package, *minba* is implemented for MaxEnt models and uses the Boyce Index as a measure of model performance. This might be extended to other algorithms and modelling techniques, as well as to other evaluation metrics.



To install the latest version:

```
library(devtools)
install_github("xavi-rp/MinBAR")
```



## References



