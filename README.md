# MinBAR

MinBAR is an R package that aims at (1) defining the minimum background extent necessary to fit Species Distribution Models reliable enough to extract ecologically relevant conclusions from them and (2) optimizing the modelling process in terms of computation demands.

The idea is to sequentially fit several concentric SDMs, with different diameter each (i.e. buffers), from the geographical centre of the species distribution to the periphery, until a satisfactory model is reached.

The main function of MinBAR is *minba()*. In the version 1.1.3 of the package, *minba()* is implemented for MaxEnt models and uses the Boyce Index as a measure of model performance. This might be extended to other algorithms and modelling techniques, as well as to other evaluation metrics.

To install the latest version:
```
library(devtools)
install_github("xavi-rp/MinBAR")
```
&nbsp;

See [this](https://cran.r-project.org/package=MinBAR/vignettes/Example_MinBAR_Balearics.html) vignette for an example of using MinBAR with a plant species from the Balearic Islands.

&nbsp;



### References

Rotllan-Puig, X. & Traveset, A. 2021. *Determining the Minimal Background Area for Species Distribution Models: MinBAR Package*. Ecological Modelling. 439:109353. \doi{10.1016/j.ecolmodel.2020.109353}
