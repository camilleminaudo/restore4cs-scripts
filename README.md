# restore4cs-scripts
This repository contains the scripts to harmonize and process field data from the RESTORE4CS sampling campaigns.
User need to specify a path to a file or folder needed to be processed, and the scripts harmonize and structure the data in a systematic manner.
Then, fluxes are being processed with different approaches inspired from the R package GoFluxYourself (https://github.com/Qepanna/GoFluxYourself/tree/master).
Data and fluxes estimates are saved locally, not on the remote.

<img src="https://github.com/camilleminaudo/restore4cs-scripts/blob/main/RESTORE4Cs_LOGO_DEF.jpg" width=50% height=50%>


## Install
First, you'll need to install the GoFluxYourself package:

To install a package from GitHub, one must first install the package
`devtools` from the CRAN:

``` r
if (!require("devtools", quietly = TRUE))
install.packages("devtools")
```

Then, install the `goFlux` package from GitHub:

``` r
try(detach("package:goFlux", unload = TRUE), silent = TRUE)
devtools::install_github("Qepanna/goFlux")
```


Then, clone the restore4cs-scripts repository, and it should be enough to have the scripts working.
Other scripts are functions needed for the code to work.


The functioning of the package depends on many other packages
(`data.table`, `dplyr`, `ggnewscale`, `ggplot2`, `graphics`,
`grDevices`, `grid`, `gridExtra`, `lubridate`, `minpack.lm`, `msm`,
`pbapply`, `plyr`, `purrr`, `rjson`, `rlist`, `SimDesign`, `stats`,
`tibble`, `tidyr`, `utils`), which will be installed when installing
`GoFluxYourself`. If prompted, it is recommended to update any
pre-installed packages.
