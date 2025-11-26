# GregoryLoredo

[![Build Status](https://github.com/stefanocovino/GregoryLoredo.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/stefanocovino/GregoryLoredo.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://stefanocovino.github.io/GregoryLoredo.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://stefanocovino.github.io/GregoryLoredo.jl/dev/)



Phil Gregory and Tom Loredo published in 1992 a [seminal paper](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..146G/abstract) about a Bayesian method to detect periodic signal in astronomical datasets without assuming a particular shape for the periodic pattern.

The alghorithm has beed widely applied in astronomical researche and this is its implementation in the `Julia` language.

An example of application is based on data published by [Raywade et al.(2020)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.495.3551R/abstract). The Authors collected [fast-radio bursts (FRBs)](https://en.wikipedia.org/wiki/Fast_radio_burst) from the event **FRB 121102** in order to identify a possible periodicity in the occurrence these events.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/stefanocovino/GregoryLoredo.jl.git")
```

will install this package.

[Here](https://stefanocovino.github.io/GregoryLoredo.jl/stable/)'s the documentation!


