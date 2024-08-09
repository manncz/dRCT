# Design-based RCT Analysis Package (dRCT)

## Installation
To install, use either the `remotes` or `devtools` package in `R` 
If `remotes` is not already installed, run:
```
> install.packages('remotes')
```

To install `dRCT`, run the following:

```
> remotes::install_github('manncz/dRCT',dep=TRUE)
```

## Usage
The main functions for estimating average effects are `loop()`, for Bernoulli experiments, and `p_loop()` for paired experiments

For detailed examples, see the tutorial materials at https://github.com/manncz/edm-rct-tutorial
