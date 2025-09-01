// Core implementation of 1-D Total Variation Denoising using Condat's 2013
// and 2017 algorithms. The routines operate directly on raw pointers and are
// exposed to Python through pybind11 with NumPy array support.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <algorithm>
#include <vector>

namespace py = pybind11;

namespace tvd {

// -----------------------------------------------------------------------------
// 2013 algorithm (version 1)
// -----------------------------------------------------------------------------

template <typename T>
void tv1d_denoise(const T *input, T *output, unsigned int width, T lambda)
{
    std::vector<unsigned int> indstart_low(width);
    std::vector<unsigned int> indstart_up(width);

    unsigned int j_low = 0, j_up = 0, jseg = 0, indjseg = 0, i = 1, indjseg2, ind;
    T output_low_first = input[0] - lambda;
    T output_low_curr = output_low_first;
    T output_up_first = input[0] + lambda;
    T output_up_curr = output_up_first;
    T twolambda = static_cast<T>(2) * lambda;

    if (width == 1)
    {
        output[0] = input[0];
        return;
    }

    indstart_low[0] = 0;
    indstart_up[0] = 0;
    width--;
    for (; i < width; ++i)
    {
        if (input[i] >= output_low_curr)
        {
            if (input[i] <= output_up_curr)
            {
                output_up_curr +=
                    (input[i] - output_up_curr) / (i - indstart_up[j_up] + 1);
                output[indjseg] = output_up_first;
                while ((j_up > jseg) &&
                       (output_up_curr <= output[ind = indstart_up[j_up - 1]]))
                    output_up_curr +=
                        (output[ind] - output_up_curr) *
                        static_cast<T>(indstart_up[j_up--] - ind) /
                        (i - ind + 1);
                if (j_up == jseg)
                {
                    while ((output_up_curr <= output_low_first) &&
                           (jseg < j_low))
                    {
                        indjseg2 = indstart_low[++jseg];
                        output_up_curr +=
                            (output_up_curr - output_low_first) *
                            static_cast<T>(indjseg2 - indjseg) /
                            (i - indjseg2 + 1);
                        while (indjseg < indjseg2)
                            output[indjseg++] = output_low_first;
                        output_low_first = output[indjseg];
                    }
                    output_up_first = output_up_curr;
                    indstart_up[j_up = jseg] = indjseg;
                }
                else
                {
                    output[indstart_up[j_up]] = output_up_curr;
                }
            }
            else
            {
                output_up_curr = output[i] = input[indstart_up[++j_up] = i];
                output_low_curr +=
                    (input[i] - output_low_curr) /
                    (i - indstart_low[j_low] + 1);
                output[indjseg] = output_low_first;
                while ((j_low > jseg) &&
                       (output_low_curr >=
                        output[ind = indstart_low[j_low - 1]]))
                    output_low_curr +=
                        (output[ind] - output_low_curr) *
                        static_cast<T>(indstart_low[j_low--] - ind) /
                        (i - ind + 1);
                if (j_low == jseg)
                {
                    while ((output_low_curr >= output_up_first) &&
                           (jseg < j_up))
                    {
                        indjseg2 = indstart_up[++jseg];
                        output_low_curr +=
                            (output_low_curr - output_up_first) *
                            static_cast<T>(indjseg2 - indjseg) /
                            (i - indjseg2 + 1);
                        while (indjseg < indjseg2)
                            output[indjseg++] = output_up_first;
                        output_up_first = output[indjseg];
                    }
                    if ((indstart_low[j_low = jseg] = indjseg) == i)
                        output_low_first = output_up_first - twolambda;
                    else
                        output_low_first = output_low_curr;
                }
                else
                {
                    output[indstart_low[j_low]] = output_low_curr;
                }
            }
        }
        else
        {
            output_up_curr +=
                ((output_low_curr = output[i] =
                      input[indstart_low[++j_low] = i]) -
                 output_up_curr) /
                (i - indstart_up[j_up] + 1);
            output[indjseg] = output_up_first;
            while ((j_up > jseg) &&
                   (output_up_curr <= output[ind = indstart_up[j_up - 1]]))
                output_up_curr +=
                    (output[ind] - output_up_curr) *
                    static_cast<T>(indstart_up[j_up--] - ind) /
                    (i - ind + 1);
            if (j_up == jseg)
            {
                while ((output_up_curr <= output_low_first) &&
                       (jseg < j_low))
                {
                    indjseg2 = indstart_low[++jseg];
                    output_up_curr +=
                        (output_up_curr - output_low_first) *
                        static_cast<T>(indjseg2 - indjseg) /
                        (i - indjseg2 + 1);
                    while (indjseg < indjseg2)
                        output[indjseg++] = output_low_first;
                    output_low_first = output[indjseg];
                }
                if ((indstart_up[j_up = jseg] = indjseg) == i)
                    output_up_first = output_low_first + twolambda;
                else
                    output_up_first = output_up_curr;
            }
            else
            {
                output[indstart_up[j_up]] = output_up_curr;
            }
        }
    }

    if (input[i] + lambda <= output_low_curr)
    {
        while (jseg < j_low)
        {
            indjseg2 = indstart_low[++jseg];
            while (indjseg < indjseg2)
                output[indjseg++] = output_low_first;
            output_low_first = output[indjseg];
        }
        while (indjseg < i)
            output[indjseg++] = output_low_first;
        output[indjseg] = input[i] + lambda;
    }
    else if (input[i] - lambda >= output_up_curr)
    {
        while (jseg < j_up)
        {
            indjseg2 = indstart_up[++jseg];
            while (indjseg < indjseg2)
                output[indjseg++] = output_up_first;
            output_up_first = output[indjseg];
        }
        while (indjseg < i)
            output[indjseg++] = output_up_first;
        output[indjseg] = input[i] - lambda;
    }
    else
    {
        output_low_curr += (input[i] + lambda - output_low_curr) /
                           (i - indstart_low[j_low] + 1);
        output[indjseg] = output_low_first;
        while ((j_low > jseg) &&
               (output_low_curr >= output[ind = indstart_low[j_low - 1]]))
            output_low_curr +=
                (output[ind] - output_low_curr) *
                static_cast<T>(indstart_low[j_low--] - ind) /
                (i - ind + 1);
        if (j_low == jseg)
        {
            if (output_up_first >= output_low_curr)
            {
                while (indjseg <= i)
                    output[indjseg++] = output_low_curr;
            }
            else
            {
                output_up_curr += (input[i] - lambda - output_up_curr) /
                                   (i - indstart_up[j_up] + 1);
                output[indjseg] = output_up_first;
                while ((j_up > jseg) &&
                       (output_up_curr <=
                        output[ind = indstart_up[j_up - 1]]))
                    output_up_curr +=
                        (output[ind] - output_up_curr) *
                        static_cast<T>(indstart_up[j_up--] - ind) /
                        (i - ind + 1);
                while (jseg < j_up)
                {
                    indjseg2 = indstart_up[++jseg];
                    while (indjseg < indjseg2)
                        output[indjseg++] = output_up_first;
                    output_up_first = output[indjseg];
                }
                indjseg = indstart_up[j_up];
                while (indjseg <= i)
                    output[indjseg++] = output_up_curr;
            }
        }
        else
        {
            while (jseg < j_low)
            {
                indjseg2 = indstart_low[++jseg];
                while (indjseg < indjseg2)
                    output[indjseg++] = output_low_first;
                output_low_first = output[indjseg];
            }
            indjseg = indstart_low[j_low];
            while (indjseg <= i)
                output[indjseg++] = output_low_curr;
        }
    }
}

// -----------------------------------------------------------------------------
// 2017 algorithm (version 2) by Laurent Condat under the CeCILL licence.
// -----------------------------------------------------------------------------

/*
Total variation denoising of 1-D signals, a.k.a. Fused lasso
signal approximator, by Laurent Condat.

Version 2.0, Aug. 30, 2017.

Usage rights : Copyright Laurent Condat.
This file is distributed under the terms of the CeCILL
licence (compatible with the GNU GPL), which can be
found at the URL "http://www.cecill.info".
*/
template <typename T>
void tv1d_denoise_v2(const T *input, T *output, unsigned int width, T lambda)
{
    std::vector<unsigned int> indstart_low(width);
    std::vector<unsigned int> indstart_up(width);

    unsigned int j_low = 0, j_up = 0, jseg = 0, indjseg = 0, i = 1, indjseg2, ind;
    T output_low_first = input[0] - lambda;
    T output_low_curr = output_low_first;
    T output_up_first = input[0] + lambda;
    T output_up_curr = output_up_first;
    T twolambda = static_cast<T>(2) * lambda;

    if (width == 1)
    {
        output[0] = input[0];
        return;
    }

    indstart_low[0] = 0;
    indstart_up[0] = 0;
    width--;
    for (; i < width; ++i)
    {
        if (input[i] >= output_low_curr)
        {
            if (input[i] <= output_up_curr)
            {
                output_up_curr +=
                    (input[i] - output_up_curr) / (i - indstart_up[j_up] + 1);
                output[indjseg] = output_up_first;
                while ((j_up > jseg) &&
                       (output_up_curr <= output[ind = indstart_up[j_up - 1]]))
                    output_up_curr +=
                        (output[ind] - output_up_curr) *
                        static_cast<T>(indstart_up[j_up--] - ind) /
                        (i - ind + 1);
                if (j_up == jseg)
                {
                    while ((output_up_curr <= output_low_first) &&
                           (jseg < j_low))
                    {
                        indjseg2 = indstart_low[++jseg];
                        output_up_curr +=
                            (output_up_curr - output_low_first) *
                            static_cast<T>(indjseg2 - indjseg) /
                            (i - indjseg2 + 1);
                        while (indjseg < indjseg2)
                            output[indjseg++] = output_low_first;
                        output_low_first = output[indjseg];
                    }
                    output_up_first = output_up_curr;
                    indstart_up[j_up = jseg] = indjseg;
                }
                else
                {
                    output[indstart_up[j_up]] = output_up_curr;
                }
            }
            else
            {
                output_up_curr = output[i] = input[indstart_up[++j_up] = i];
                output_low_curr +=
                    (input[i] - output_low_curr) /
                    (i - indstart_low[j_low] + 1);
                output[indjseg] = output_low_first;
                while ((j_low > jseg) &&
                       (output_low_curr >=
                        output[ind = indstart_low[j_low - 1]]))
                    output_low_curr +=
                        (output[ind] - output_low_curr) *
                        static_cast<T>(indstart_low[j_low--] - ind) /
                        (i - ind + 1);
                if (j_low == jseg)
                {
                    while ((output_low_curr >= output_up_first) &&
                           (jseg < j_up))
                    {
                        indjseg2 = indstart_up[++jseg];
                        output_low_curr +=
                            (output_low_curr - output_up_first) *
                            static_cast<T>(indjseg2 - indjseg) /
                            (i - indjseg2 + 1);
                        while (indjseg < indjseg2)
                            output[indjseg++] = output_up_first;
                        output_up_first = output[indjseg];
                    }
                    if ((indstart_low[j_low = jseg] = indjseg) == i)
                        output_low_first = output_up_first - twolambda;
                    else
                        output_low_first = output_low_curr;
                }
                else
                {
                    output[indstart_low[j_low]] = output_low_curr;
                }
            }
        }
        else
        {
            output_up_curr +=
                ((output_low_curr = output[i] =
                      input[indstart_low[++j_low] = i]) -
                 output_up_curr) /
                (i - indstart_up[j_up] + 1);
            output[indjseg] = output_up_first;
            while ((j_up > jseg) &&
                   (output_up_curr <= output[ind = indstart_up[j_up - 1]]))
                output_up_curr +=
                    (output[ind] - output_up_curr) *
                    static_cast<T>(indstart_up[j_up--] - ind) /
                    (i - ind + 1);
            if (j_up == jseg)
            {
                while ((output_up_curr <= output_low_first) &&
                       (jseg < j_low))
                {
                    indjseg2 = indstart_low[++jseg];
                    output_up_curr +=
                        (output_up_curr - output_low_first) *
                        static_cast<T>(indjseg2 - indjseg) /
                        (i - indjseg2 + 1);
                    while (indjseg < indjseg2)
                        output[indjseg++] = output_low_first;
                    output_low_first = output[indjseg];
                }
                if ((indstart_up[j_up = jseg] = indjseg) == i)
                    output_up_first = output_low_first + twolambda;
                else
                    output_up_first = output_up_curr;
            }
            else
            {
                output[indstart_up[j_up]] = output_up_curr;
            }
        }
    }

    if (input[i] + lambda <= output_low_curr)
    {
        while (jseg < j_low)
        {
            indjseg2 = indstart_low[++jseg];
            while (indjseg < indjseg2)
                output[indjseg++] = output_low_first;
            output_low_first = output[indjseg];
        }
        while (indjseg < i)
            output[indjseg++] = output_low_first;
        output[indjseg] = input[i] + lambda;
    }
    else if (input[i] - lambda >= output_up_curr)
    {
        while (jseg < j_up)
        {
            indjseg2 = indstart_up[++jseg];
            while (indjseg < indjseg2)
                output[indjseg++] = output_up_first;
            output_up_first = output[indjseg];
        }
        while (indjseg < i)
            output[indjseg++] = output_up_first;
        output[indjseg] = input[i] - lambda;
    }
    else
    {
        output_low_curr += (input[i] + lambda - output_low_curr) /
                           (i - indstart_low[j_low] + 1);
        output[indjseg] = output_low_first;
        while ((j_low > jseg) &&
               (output_low_curr >= output[ind = indstart_low[j_low - 1]]))
            output_low_curr +=
                (output[ind] - output_low_curr) *
                static_cast<T>(indstart_low[j_low--] - ind) /
                (i - ind + 1);
        if (j_low == jseg)
        {
            if (output_up_first >= output_low_curr)
            {
                while (indjseg <= i)
                    output[indjseg++] = output_low_curr;
            }
            else
            {
                output_up_curr += (input[i] - lambda - output_up_curr) /
                                   (i - indstart_up[j_up] + 1);
                output[indjseg] = output_up_first;
                while ((j_up > jseg) &&
                       (output_up_curr <=
                        output[ind = indstart_up[j_up - 1]]))
                    output_up_curr +=
                        (output[ind] - output_up_curr) *
                        static_cast<T>(indstart_up[j_up--] - ind) /
                        (i - ind + 1);
                while (jseg < j_up)
                {
                    indjseg2 = indstart_up[++jseg];
                    while (indjseg < indjseg2)
                        output[indjseg++] = output_up_first;
                    output_up_first = output[indjseg];
                }
                indjseg = indstart_up[j_up];
                while (indjseg <= i)
                    output[indjseg++] = output_up_curr;
            }
        }
        else
        {
            while (jseg < j_low)
            {
                indjseg2 = indstart_low[++jseg];
                while (indjseg < indjseg2)
                    output[indjseg++] = output_low_first;
                output_low_first = output[indjseg];
            }
            indjseg = indstart_low[j_low];
            while (indjseg <= i)
                output[indjseg++] = output_low_curr;
        }
    }
}

// -----------------------------------------------------------------------------
// Wrappers operating on NumPy arrays
// -----------------------------------------------------------------------------

template <typename T>
py::array_t<T> TVD(py::array_t<T, py::array::c_style | py::array::forcecast> in,
                   double lambda)
{
    auto buf = in.request();
    auto n = static_cast<unsigned int>(buf.size);
    const T *data = static_cast<T *>(buf.ptr);
    py::array_t<T> result(buf.size);
    tv1d_denoise(data, result.mutable_data(), n, static_cast<T>(lambda));
    return result;
}

template <typename T>
py::array_t<T> TVD_v2(py::array_t<T, py::array::c_style | py::array::forcecast> in,
                      double lambda)
{
    auto buf = in.request();
    auto n = static_cast<unsigned int>(buf.size);
    const T *data = static_cast<T *>(buf.ptr);
    py::array_t<T> result(buf.size);
    tv1d_denoise_v2(data, result.mutable_data(), n, static_cast<T>(lambda));
    return result;
}

template <typename T>
py::array_t<T> D_TVD_R(py::array_t<T, py::array::c_style | py::array::forcecast> in,
                        double lambda)
{
    auto buf = in.request();
    auto n = static_cast<unsigned int>(buf.size);
    const T *data = static_cast<T *>(buf.ptr);

    std::vector<T> detrend(n), denoised(n);

    T base_value = data[0];
    T mean = 0;
    size_t sample = static_cast<size_t>(n / 100);
    if (sample < 1)
        sample = 1;
    else if (sample > 10)
        sample = 10;
    for (size_t k = 0; k < sample; ++k)
        mean += data[k];
    mean /= static_cast<T>(sample);
    base_value = mean;

    detrend[0] = base_value;
    for (unsigned int idx = 1; idx < n; ++idx)
        detrend[idx] = data[idx] - data[idx - 1];

    tv1d_denoise(detrend.data(), denoised.data(), n, static_cast<T>(lambda));

    py::array_t<T> result(n);
    T *out = result.mutable_data();
    for (unsigned int idx = 0; idx < n; ++idx)
    {
        if (idx == 0)
            out[idx] = base_value;
        else
            out[idx] = denoised[idx] + out[idx - 1];
    }
    return result;
}

template <typename T>
py::array_t<T> D_TVD_R_v2(py::array_t<T, py::array::c_style | py::array::forcecast> in,
                           double lambda)
{
    auto buf = in.request();
    auto n = static_cast<unsigned int>(buf.size);
    const T *data = static_cast<T *>(buf.ptr);

    std::vector<T> detrend(n), denoised(n);

    T base_value = data[0];
    T mean = 0;
    size_t sample = static_cast<size_t>(n / 100);
    if (sample < 1)
        sample = 1;
    else if (sample > 10)
        sample = 10;
    for (size_t k = 0; k < sample; ++k)
        mean += data[k];
    mean /= static_cast<T>(sample);
    base_value = mean;

    detrend[0] = base_value;
    for (unsigned int idx = 1; idx < n; ++idx)
        detrend[idx] = data[idx] - data[idx - 1];

    tv1d_denoise_v2(detrend.data(), denoised.data(), n,
                     static_cast<T>(lambda));

    py::array_t<T> result(n);
    T *out = result.mutable_data();
    for (unsigned int idx = 0; idx < n; ++idx)
    {
        if (idx == 0)
            out[idx] = base_value;
        else
            out[idx] = denoised[idx] + out[idx - 1];
    }
    return result;
}

} // namespace tvd

PYBIND11_MODULE(TVDCondat2013, m)
{
    m.doc() = R"pbdoc(
        TVDCondat2013 provides 1-D Total Variation Denoising routines
        based on Condat's algorithms.
    )pbdoc";

    m.def("TVD", &tvd::TVD<double>, py::arg("signal").noconvert(),
          py::arg("lambda"),
          "Total variation denoising for float64 arrays");
    m.def("TVD", &tvd::TVD<float>, py::arg("signal").noconvert(),
          py::arg("lambda"),
          "Total variation denoising for float32 arrays");

    m.def("D_TVD_R", &tvd::D_TVD_R<double>, py::arg("signal").noconvert(),
          py::arg("lambda"),
          "Detrend-denoise-retrend for float64 arrays");
    m.def("D_TVD_R", &tvd::D_TVD_R<float>, py::arg("signal").noconvert(),
          py::arg("lambda"),
          "Detrend-denoise-retrend for float32 arrays");

    m.def("TVD_v2", &tvd::TVD_v2<double>, py::arg("signal").noconvert(),
          py::arg("lambda"),
          "Total variation denoising (algorithm v2) for float64 arrays");
    m.def("TVD_v2", &tvd::TVD_v2<float>, py::arg("signal").noconvert(),
          py::arg("lambda"),
          "Total variation denoising (algorithm v2) for float32 arrays");

    m.def("D_TVD_R_v2", &tvd::D_TVD_R_v2<double>,
          py::arg("signal").noconvert(), py::arg("lambda"),
          "Detrend-denoise-retrend (algorithm v2) for float64 arrays");
    m.def("D_TVD_R_v2", &tvd::D_TVD_R_v2<float>,
          py::arg("signal").noconvert(), py::arg("lambda"),
          "Detrend-denoise-retrend (algorithm v2) for float32 arrays");
}

