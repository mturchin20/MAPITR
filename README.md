# MAPIT_R: MArginal EPIstasis Test for Regions

[![CRAN status badge](https://www.r-pkg.org/badges/version/MAPIT_R)](https://cran.r-project.org/package=MAPIT_R)
[![Travis Build Status](https://travis-ci.org/mturchin20/MAPIT_R.svg?branch=master)](https://travis-ci.org/mturchin20/MAPIT_R)

The `MAPIT_R` R package provides accessible functions for running the
algorithms described in...

If you find a bug, or you have a question or feedback on our work,
please post an [issue][issues].

## Citing this work

If you find the `MAPIT_R` package or any of the source code in this
repository useful for your work, please cite:

> Turchin MC, Darnell G, Crawford L, and Ramachandran S (2020) 
> "...

## License

Copyright (c) 2020, Michael Turchin, Sohini Ramachandran, and Lorin Crawford.

All source code and software in this repository are made available
under the terms of the [MIT license][mit-license]. See
file [LICENSE](LICENSE) for the full text of the license.

## Quick Start

To install `MAPIT_R` from [CRAN](https://cran.r-project.org/web/packages/MAPIT_R/index.html):

```{r}
install.packages("MAPIT_R")
```

To install the most recent dev version of `MAPIT_R` from [github](https://github.com/mturchin20/MAPIT_R):
```{r}
install.packages("devtools")
devtools::install_github("mturchin20/MAPIT_R@v1.0.3", build_vignettes=TRUE)
```

Once you have installed the package, load the package in R:

```{r}
library("MAPIT_R")
```

Next, view and run the example code provided in the 
[introductory vignette][MAPIT_R-vignette1] using simulated data. 

## Credits

The `MAPIT_R` R package was developed by [Michael Turchin][michaelt] at
the [Brown University][brownu], with contributions from
[Greg Darnell][greg] and [Lorin Crawford][lorin].

[MAPIT_R-website]: http://mturchin20.github.io/MAPIT_R/ 
[MAPIT_R-vignette1]: http://mturchin20.github.io/MAPIT_R/articles/MAPIT_RIntro.1.SimulatedData.html
[MAPIT_R-vignette2]: http://mturchin20.github.io/MAPIT_R/articles/MAPIT_RIntro.2.RealData.html
[biorxiv-paper]: 
[issues]: https://github.com/mturchin20/MAPIT_R/issues
[lorin]: http://www.lcrawlab.com/ 
[michaelt]: http://home.uchicago.edu/mturchin20/index.html 
[mit-license]: https://opensource.org/licenses/mit-license.html
[greg]: https://www.gregdarnell.com/ 
[uchicago]: https://www.brown.edu
