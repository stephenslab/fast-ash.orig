# fast-ash

The goal of this repo is to look at methods for making ashr run faster, possibly at cost of accuracy. 
Motivation is particularly for 
large data sets and when ash is used repeatedly as part of an iterative method (so each application of ash does not
necessarily need to run to convergence).

Some things that affect speed are:
- number of samples (n)
- grid size (k)
- convergence tolerance
- starting point? (if EM used)

Initial ideas include to drop grid points when they don't seem to get any weight, and analyze only subsets of data. Maybe in combination with some initial low-tolerance to reduce time.

# Contributing 

This repo follows the [ashlar](http://github.com/stephenslab/ashlar) template.

Add analyses in the 'analysis' subdirectory.








