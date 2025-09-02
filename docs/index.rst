TVDCondat2013 Documentation
=======================================

TVDCondat2013 provides fast 1‑D total variation denoising routines based on
Laurent Condat's algorithms. The core is written in C++ and exposed to Python
through pybind11, offering a tiny, dependency‑free extension.  Two
implementations are available:

* ``TVD`` and ``D_TVD_R`` — direct ports of the 2013 algorithm.
* ``TVD_v2`` and ``D_TVD_R_v2`` — an accelerated variant derived from Condat's
  2017 work with improved convergence on large signals.

Quick start
-----------

.. code-block:: python

    import numpy as np
    from TVDCondat2013 import TVD, TVD_v2

    noisy = np.random.randn(100)
    denoised = TVD(noisy, 10.0)       # 2013 algorithm
    denoised_v2 = TVD_v2(noisy, 10.0) # 2017 algorithm

``D_TVD_R`` and ``D_TVD_R_v2`` offer the same pair of algorithms for the
detrend–denoise–retrend workflow.  See the :doc:`api` reference for the complete
list of functions and their arguments.

References
----------

* L. Condat, "A Direct Algorithm for 1D Total Variation Denoising," *IEEE
  Signal Processing Letters*, 2013.
* L. Condat, "Fast Projection onto the Simplex and the L1 Ball," *Mathematical
  Programming*, 2017.

.. toctree::
   :maxdepth: 1

   api

