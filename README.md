## What is `climate4R.climdex`?

**climate4R** is a set of R packages for transparent climate data access, post processing (including bias correction and downscaling) and visualization. For more information and references, visit the [climate4R page](http://www.meteo.unican.es/climate4r).

`climate4R.climdex` is a wrapper of the R package [`climdex.pcic`](https://github.com/pacificclimate/climdex.pcic) allowing the calculation of the [ETCCDI core indices](http://etccdi.pacificclimate.org/list_27_indices.shtml) for a seamless integration with the **climate4R** data structures, and providing support for parallel computing.


****

### Installation

The recommended procedure for installing the package is using the devtools package. Note that this package depends on [`transformeR`](https://github.com/SantanderMetGroup/transformeR), another package from the **climate4R** bundle. Thus:

```R
devtools::install_github(c("SantanderMetGroup/transformeR", "SantanderMetGroup/climate4R.climdex"))
```

**** 

### To start...

```R
library(climate4R.climdex)
```

An overview of the `climdex.pcic`package functionalities can be obtained using

```R
?climdex.pcic
```

A list of all available indices, the atomic functions calculating them, the input variables required... is printed on screen with:

```R
climdexShow()
```

See the examples in: 

```R
?climdexGrid
``` 

