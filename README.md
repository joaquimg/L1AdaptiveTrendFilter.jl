# L1AdaptiveTrendFilter.jl

L1AdaptiveTrendFilter.jl is a Julia package for filtering and identifying unknown underlying trends of a given noisy signal. It is originally a supplemental material for our paper ["L1 Adaptive Trend Filter via Fast Coordinate Descent"][l1adaptivepaper].

## Installing

In order to use it, download the package in the julia prompt with the command:
```
  Pkg.clone("git://github.com/joaquimg/L1AdaptiveTrendFilter.jl.git")
```

Then simply include the command `using L1AdaptiveTrendFilter` to import the package.

## Filtering

So far there is only one function: `l1_adaptive_trend_filter` 

The inputs are described below:

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
