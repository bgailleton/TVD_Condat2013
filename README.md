TVDCondat2013
==============

*UPDATE 2025* : Still using that code as is, I updated the references and I am now preparing a `pip` package precompiled for Window/Mac/Linux. I'm also updating the method to a more optimised version.

Version 2 of the tools implements Condat's 2017 improvements, providing faster convergence while keeping the same interface.

TVDCondat2013 is a python portage of the 1D Total Variation Denoising algorithm from Condat 2013: _A Direct Algorithm for 1D Total Variation Denoising_ (Sign. Proc. Letters, DOI:10.1109/LSP.2013.2278339) using pybind11 to bind C++ and NumPy directly.

The `c++` core code has been adapted from `C` version available on the website of the manuscript orignal authors: [The Paper](https://lcondat.github.io/publis/Condat-fast_TV-SPL-2013.pdf)
*Cite it if you use it.*

This package is mostly to train myself packging a python package from `c++` but also it is a really useful and efficient algorithm for:
- Direct denoising for data that can be represented by flat segments
- Indirect curve denoising using a detrend-denoise-retrend approach (Not implemented yet)

*The package still is under active development*

This work

Quick start
------------

So far the following denoising of a `numpy` array are implemented. The
original 2013 algorithm, the faster 2017 variant, a taut string algorithm, and
a fused lasso variant with an additional :math:`\ell_1` penalty are exposed:

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

More working examples in the `examples` folder.

Installation
------------

I am working on developing a `pip` and a `conda` package at point.

**On Unix (Linux, OS X)**

 - clone this repository
 - `pip install ./TVDCondat2013`

**On Windows (Requires Visual Studio 2015)**

 - For Python 3.5:
     - clone this repository
     - `pip install ./TVDCondat2013`
 - For earlier versions of Python, including Python 2.7:

   Building requires a C++14 compliant compiler (i.e. Visual Studio 2015 on
   Windows). Running a regular `pip install` command will detect the version
   of the compiler used to build Python and attempt to build the extension
   with it. We must force the use of Visual Studio 2015.

     - clone this repository
     - `"%VS140COMNTOOLS%\..\..\VC\vcvarsall.bat" x64`
     - `set DISTUTILS_USE_SDK=1`
     - `set MSSdk=1`
     - `pip install ./TVDCondat2013`

   Note that this requires the user building `TVDCondat2013` to have registry edition
   rights on the machine, to be able to run the `vcvarsall.bat` script.


Windows runtime requirements
----------------------------

On Windows, the Visual C++ 2015 redistributable packages are a runtime
requirement for this project. It can be found [here](https://www.microsoft.com/en-us/download/details.aspx?id=48145).

If you use the Anaconda python distribution, you may require the Visual Studio
runtime as a platform-dependent runtime requirement for you package:

```yaml
requirements:
  build:
    - python
    - setuptools
    - pybind11

  run:
   - python
   - vs2015_runtime  # [win]
```


Building the documentation
--------------------------

Documentation for the example project is generated using Sphinx. Sphinx has the
ability to automatically inspect the signatures and documentation strings in
the extension module to generate beautiful documentation in a variety formats.
The following command generates HTML-based reference documentation; for other
formats please refer to the Sphinx manual:

 - `TVDCondat2013/docs`
 - `pip install sphinx furo`
 - `make html`


Running the tests
-----------------

Running the tests requires `pytest`.

```bash
py.test .
```

Wheel integrity checks
----------------------

The repository ships with `wheelhouse/debug_wheel.py` which inspects wheel
archives for CRC mismatches and shifted local headers. If the script reports a
corrupted member it means the payload was modified after the wheel was built
(`strip` and some antivirus scanners are known to do this). Always rebuild the
wheel instead of editing an existing archive in place; otherwise `pip` will
reject the file. The release workflow runs this validator on the freshly
produced wheels and fails fast if corruption is detected. Because every run
builds wheels from a clean checkout, removing any committed artefacts and
rerunning the workflow is sufficient to resolve earlier corruption reports.

Copyright
---------

Portions of this software are Copyright (c) 2013-2025 Laurent Condat.
The taut string algorithm is adapted from Matlab code by Lutz Dümbgen.
