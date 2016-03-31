# L1AdaptiveTrendFilter.jl

L1AdaptiveTrendFilter.jl is a Julia package for filtering and identifying unknown underlying trends of a given noisy signal. It is originally a supplemental material for our paper ["L1 Adaptive Trend Filter via Fast Coordinate Descent"][l1adaptivepaper].

## Installing

In order to use it, download the package in the julia prompt with the command:
```
Pkg.clone("git://github.com/joaquimg/L1AdaptiveTrendFilter.jl.git")
```

Then simply include the command `using L1AdaptiveTrendFilter` to import the package.

## Filtering

The function `l1_adaptive_trend_filter` performs the filtering via coordinate descent and takes the following required inputs:

* y : Signal or time-series to be filtered;
* components: List of integers corresponding to the types of components to be considered.

And some optional inputs:

* f: Vector of overcomplete frequencies (Default=ø);
* numλ: Size of the regularizer path for the parameter λ (Default=40);
* numγ: Size of the regularizer path for the parameter γ (Default=10);
* MAXITER: Maximum number of iterations (Default=500);
* verbose: Boolean flag for displaying progress of algorithm (Default=false);
* lower_bounds: List of lower bounds for each component type (Default=[-∞,-∞,-∞,-∞,-∞]);
* upper_bounds: List of upper bpunds for each component type (Default=[+∞,+∞,+∞,+∞,+∞]).

This function returns:

*  β_path: Path of components coefficients;
*  y_path: Path of filtered signals;
*  β_best: Best components coefficients according to the EBIC criteria;
*  y_best: Best filtered signal according to the EBIC criteria;
*  λ_best: Best value for the λ regularizer according to the EBIC criteria;
*  γ_best: Best value for the γ regularizer according to the EBIC criteria.

## Example

## Citing this package

Please use the following Bibtex reference:
```
@article{souto2016ell_1,
  title={$$\backslash$ ell\_1 $ Adaptive Trend Filter via Fast Coordinate Descent},
  author={Souto, Mario and Garcia, Joaquim D and Amaral, Gustavo C},
  journal={arXiv preprint arXiv:1603.03799},
  year={2016}
}
```

[l1adaptivepaper]: http://arxiv.org/pdf/1603.03799.pdf
