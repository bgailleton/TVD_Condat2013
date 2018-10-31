TVDCondat2013
==============

TVDCondat2013 is a python portage of the 1D Total Variation Denoising algorithm from Condat 2013: _A Direct Algorithm for 1D Total Variation Denoising_ (Sign. Proc. Letters, DOI:10.1109/LSP.2013.2278339) using xtensor and pybind11 to bind c++ and numpy. 

The `c++` core code has been adapted from `C` version available on the website of the manuscript orignal authors: http://www.gipsa-lab.grenoble-inp.fr/~laurent.condat/publications.html
*Cite it if you use it.*

This package is mostly to train myself packging a python package from `c++` but also it is a really useful and efficient algorithm for:
- Direct denoising for data that can be represented by flat segments
- Indirect curve denoising using a detrend-denoise-retrend approach (Not implemented yet)

*The package still is under active development*

Quick start
------------

So far the denoising of a `numpy` array is implemented:

*Quick use*
```
from TVDCondat2013 import TVD
...
denoised = TVD(MyNumpyArray,lambda_TVD)
...
```

![Effect of regulation parameter lambda on the TVD](https://raw.githubusercontent.com/bgailleton/TVD_Condat2013/master/examples/Example.png)

*Full example with plotting:*
```
import numpy as np
from TVDCondat2013 import TVD
from matplotlib import pyplot as plt


# Generating n segments of a noisy signal
size_segment = 400
n_segments = 5
C = np.random.rand(size_segment) + np.random.randint(-5,5)
for d in range(n_segments-1):
  A = np.random.rand(size_segment) + np.random.randint(-5,5)
  C = np.concatenate((C,A))

X = np.arange(n_segments*size_segment)

plt.plot(X,C,color= 'r', lw = 0.2, zorder = 1, label= "Raw signal")

# Setting the regulation parameters
lambda_TVD = [0.1,1,10,100,1000]

# Denoising with different lambda

for l in lambda_TVD:
  denoised = TVD(C,l)
  plt.plot(X,denoised, lw = 1, zorder = 2, label= r"TVD $\lambda$ = %s"%(l))



plt.xlabel("X")
plt.ylabel("Signal")

plt.legend()

plt.savefig("Example.png", dpi = 500)


```


Installation
------------

This guide is directly from `xtensor` documentation, let me know if this doesn't work.
I am working on developping a `pip` and a `conda` package at point.

**On Unix (Linux, OS X)**

 - clone this repository
 - `pip install ./TVDCondat2013`

**On Windows (Requires Visual Studio 2015)**

 - For Python 3.5:
     - clone this repository
     - `pip install ./TVDCondat2013`
 - For earlier versions of Python, including Python 2.7:

   xtensor requires a C++14 compliant compiler (i.e. Visual Studio 2015 on
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
