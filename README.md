TVDCondat2013
==============

TVDCondat2013 is a python portage of the 1D Total Variation Denoising algorithm from Condat 2013 and subsequent improvements: _A Direct Algorithm for 1D Total Variation Denoising_ (Sign. Proc. Letters, DOI:10.1109/LSP.2013.2278339) using pybind11 to bind C++ and NumPy directly.

This is a minimal implementation with `pybind11`


The `c++` core code has been adapted from `C` version available on the website of the manuscript orignal authors: [The Paper](https://lcondat.github.io/publis/Condat-fast_TV-SPL-2013.pdf)


*Cite it if you use it.*



Straightfoward API
------------

The following denoising of a `numpy` array are implemented. The
original 2013 algorithm, the faster 2017 variant, a taut string algorithm, and
a fused lasso variant with an additional :math:`\ell_1` penalty are exposed:

generic api: `method(numpy_array_1D, scalar_regulator)`

```
from TVDCondat2013 import fused_lasso, tvd_2013, tvd_2017, tvd_tautstring
...
denoised_v1 = tvd_2013(MyNumpyArray, lambda_TVD)       # 2013 algorithm
denoised_v2 = tvd_2017(MyNumpyArray, lambda_TVD)       # 2017 algorithm
denoised_ts = tvd_tautstring(MyNumpyArray, lambda_TVD) # taut string
denoised_fl = fused_lasso(MyNumpyArray, lambda_TVD, mu_L1) # fused lasso
...

```

Run `python examples/example_readme.py` to generate a figure comparing the
original, noisy, and denoised signals using all algorithms. The script builds a
piecewise-constant signal, corrupts it with Gaussian noise, denoises it with
``tvd_2013``, ``tvd_2017``, ``tvd_tautstring`` and ``fused_lasso``, and plots
the results in aligned subplots, saving the figure as `examples/Example.png`.


Installation
------------

Using `pip`:

`pip install tvdcondat2013`

Binaries are CI using github action for `python` 3.9-3.13, for Linux, Windows (64 bits) and MacOS (both ARM and intel).

From source, install `numpy, pybind11, scikit-build`, clone this repo and `pip install .` in the directory.

Copyright
---------

Portions of this software are Copyright (c) 2013-2025 Laurent Condat.
The taut string algorithm is adapted from Matlab code by Lutz Dümbgen.

Full info [here](https://lcondat.github.io/software.html)

License:
--------

CeCILL 2.1: strong copyleft open-source

Author:
-------

Boris Gailleton - boris.gailleton@univ-rennes.fr

Université de Rennes - University of Edinburgh - GFZ Potsdam
