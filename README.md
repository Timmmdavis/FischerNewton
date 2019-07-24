# FischerNewton

[![Build Status](https://travis-ci.com/Timmmdavis/FischerNewton.svg?token=1HhESyMNyqzV8R22Pqq6&branch=master)](https://travis-ci.com/Timmmdavis/FischerNewton/)

A Julia version of NUM4LCP.m (fischer_newton.m) in: 
[CutAndDisplace](https://github.com/Timmmdavis/CutAndDisplace)

This code is:
* Memory efficient (allocations in main loop are only due to: Krylov solvers & view statements )
* Relativly easy to call/adapt 
* Contains some basic regression tests (comparing directly to MATLAB implementation results). 

Original code by [Kenny Erleben](https://github.com/erleben/num4lcp/wiki/Welcome-to-the-num4lcp-wiki)

This is faster than the MATLAB version of FischerNewton in [CutAndDisplace](https://github.com/Timmmdavis/CutAndDisplace):

| NameOfFile | Size of A | MATLAB (seconds)  | Julia (seconds) | Relative speedup |
| ------------- |:----------------:| -----------------:| -------------:  | --------------:  |
| FischerNewton\test\Matricies.mat     | 1560*1560     | 1.2   |  0.8  |  1.5  |
| FischerNewton\test\Matricies_1000.mat | 5225*5225    | 22.3   | 19.5  |  1.14  |


