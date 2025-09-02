TVDCondat2013 Documentation
=======================================

TVDCondat2013 provides fast 1-D total variation denoising routines based on
Laurent Condat's algorithms. The core is written in C++ and exposed to Python
through pybind11, offering a tiny, dependency-free extension.

Quick start
-----------

.. code-block:: python

    import numpy as np
    from TVDCondat2013 import TVD

    noisy = np.random.randn(100)
    denoised = TVD(noisy, 10.0)

The module also exposes ``TVD_v2`` and ``D_TVD_R`` helpers for alternative
use-cases. See the :doc:`api` reference for the complete list of functions and
their arguments.

References
----------

* L. Condat, "A Direct Algorithm for 1D Total Variation Denoising," *IEEE
  Signal Processing Letters*, 2013.
* L. Condat, "Fast Projection onto the Simplex and the L1 Ball," *Mathematical
  Programming*, 2017.

.. toctree::
   :maxdepth: 1

   api

