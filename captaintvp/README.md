# captaintvp

Univariate trend plus autoregression routines from the CAPTAIN `univdemo`
air-passenger example.

The exported functions follow the MATLAB argument order and names:

- `fcast(y, n)`
- `aic(y, mod, out)`
- `univopt(y, p, TVP, nvr0, tar, out, Int)`
- `univ(y, ARp, TVP, nvr, ARt, Int, sm)`

Kalman filtering and fixed-interval smoothing are provided by
[KFAS](https://cran.r-project.org/package=KFAS). The airline series is
`datasets::AirPassengers`, the same 144 monthly values as CAPTAIN's `air.dat`.
IAR/DIAR trend polynomials and variance intervention are accepted by the
interface and report that they are not implemented.

See `vignette("univdemo-air", package = "captaintvp")`.
