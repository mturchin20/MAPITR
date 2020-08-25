# MAPITR: MArginal ePIstasis Test for Regions

[![CRAN status badge](https://www.r-pkg.org/badges/version/MAPITR)](https://cran.r-project.org/package=MAPITR)
[![Travis Build Status](https://travis-ci.org/mturchin20/MAPITR.svg?branch=master)](https://travis-ci.org/mturchin20/MAPITR)

The `MAPITR` R package provides accessible functions for running the
algorithms described in...

If you find a bug, or you have a question or feedback on our work,
please post an [issue][issues].

## Citing this work

If you find the `MAPITR` package or any of the source code in this
repository useful for your work, please cite:

> Turchin MC, Darnell G, Crawford L, and Ramachandran S (2020) 
> "Pathway Analysis within Multiple Human Ancestries Reveals 
> Novel Signals for Epistasis in Complex Traits".

## License

Copyright (c) 2020, Michael Turchin, Gregory Darnell, Lorin Crawford, and Sohini Ramachandran.

All source code and software in this repository are made available
under the terms of the [MIT license][mit-license]. See
file [LICENSE](LICENSE) for the full text of the license.

## Quick Start

To install `MAPITR` from [CRAN](https://cran.r-project.org/web/packages/MAPITR/index.html):

```{r}
install.packages("MAPITR")
```

To install the most recent dev version of `MAPITR` from [github](https://github.com/mturchin20/MAPITR):
```{r}
install.packages("devtools")
devtools::install_github("mturchin20/MAPITR@v1.0.4", build_vignettes=TRUE)
```

Once you have installed the package, load the package in R:

```{r}
library("MAPITR")
```

Next, view and run the example code provided in the 
[introductory vignette][MAPITR-vignette1] using simulated data. 

## Credits

The `MAPITR` R package was developed by [Michael Turchin][michaelt] at [Brown University][brownu], with contributions from
[Peter Carbonetto][peter] and [Lorin Crawford][lorin].

[MAPITR-website]: http://mturchin20.github.io/MAPITR 
[MAPITR-vignette1]: http://mturchin20.github.io/MAPITR/articles/MAPITR.Intro.SimulatedData.html
[biorxiv-paper]: http://www.google.com 
[issues]: https://github.com/mturchin20/MAPITR/issues
[lorin]: http://www.lcrawlab.com/ 
[michaelt]: http://home.uchicago.edu/mturchin20/index.html 
[mit-license]: https://opensource.org/licenses/mit-license.html
[peter]: https://pcarbo.github.io/
[uchicago]: https://www.brown.edu
