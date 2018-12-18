## What is `climate4R.climdex`?

**climate4R** is a set of R packages for transparent climate data access, post processing (including bias correction and downscaling) and visualization. For more information and references, visit the [climate4R page](http://www.meteo.unican.es/climate4r).

`climate4R.climdex` is a wrapper of the R package [`climdex.pcic`](https://github.com/pacificclimate/climdex.pcic) allowing the calculation of the [ETCCDI core indices](http://etccdi.pacificclimate.org/list_27_indices.shtml) for a seamless integration with the **climate4R** data structures, and providing support for parallel computing.


****

### Installation

The recommended procedure for installing the package is using the devtools package. Note that this package depends on [`transformeR`](https://github.com/SantanderMetGroup/transformeR), another package from the **climate4R** bundle. Thus:

```R
devtools::install_github(c("SantanderMetGroup/transformeR", "SantanderMetGroup/climate4R.climdex"))
```

A list of all available indices and the atomic functions calculating them is printed on screen with:

```R
library(climate4R.climdex)
climdexShow()
?climdexGrid   # see the examples 
```

Reference and further information: 

**[General description of the climate4R framework]** Iturbide et al. (2019) The R-based climate4R open framework for reproducible climate data access and post-processing. **Environmental Modelling and Software**, 111, 42-54. https://doi.org/10.1016/j.envsoft.2018.09.009
Check out the companion notebooks for the two examples [GitHub](https://github.com/SantanderMetGroup/notebooks).

