# CL-GAM (R package `clgam`)

Composite-link GAMMs for nested areal counts. The installable package is
[`clgam/`](clgam/) (version **0.1.25**).

## Install

```r
remotes::install_github("idaejin/CL-GAM", subdir = "clgam")
```

From a local clone of this repository:

```r
install.packages("clgam", repos = NULL, type = "source")
```

From inside `clgam/`:

```r
install.packages(".", repos = NULL, type = "source")
```

API, formula interface (`s(x1, x2, ndx, bdeg, pord)`), and examples:
[`clgam/README.md`](clgam/README.md).

## License

GPL-2
