TVDCondat2013
==============

*UPDATE 2025* : Still using that code as is, I updated the references and I am now preparing a `pip` package precompiled for Window/Mac/Linux. I'm also updating the method to a more optimised version. 

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

So far the following denoising of a `numpy` array are implemented:

**Quick use of the original denoising**
```
from TVDCondat2013 import TVD
...
denoised = TVD(MyNumpyArray,lambda_TVD)
...

```

Run `python examples/example_readme.py` to generate a figure comparing the original, noisy, and TVD-denoised signals. The script builds a piecewise-constant signal, corrupts it with Gaussian noise, denoises it with TVD, and plots the three signals in aligned subplots, saving the result as `examples/Example.png`.

**More experimental: curve denoising. So far the boundary condition might shift up or down the data. I am working on it**

```
from TVDCondat2013 import D_TVD_R
...
curve_denoised = D_TVD_R((MyNumpyArray_of_curve,lambda_TVD))
...

```

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
 - `make html`


Running the tests
-----------------

Running the tests requires `pytest`.

```bash
py.test .
```
