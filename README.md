# GregoryLoredo

[![Build Status](https://github.com/stefanocovino/GregoryLoredo.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/stefanocovino/GregoryLoredo.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://stefanocovino.github.io/GregoryLoredo.jl/stable/)



Phil Gregory and Tom Loredo published in 1992 a [seminal paper](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..146G/abstract) about a Bayesian method to detect periodic signal in astronomical datasets without assuming a particular shape for the periodic pattern.

The alghorithm has beed widely applied in astronomical research and this is its implementation in the `Julia` language. At present only the case of arrival times affedted by Poissonian statistics is implemented. The case of flux data affected by Gaussian statistics, introduced in [Gregory (1999)](https://ui.adsabs.harvard.edu/abs/1999ApJ...520..361G/abstract) will be implemented in a future release.


## Installation

The package is not registered yet, so you need the full path:

```julia
using Pkg
Pkg.add(url="https://github.com/stefanocovino/GregoryLoredo.jl.git")
```

will install this package.

[Here](https://stefanocovino.github.io/GregoryLoredo.jl/stable/)'s the documentation!


