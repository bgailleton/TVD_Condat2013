// Core implementation of 1-D Total Variation Denoising using Condat's 2013
// algorithm. The heavy lifting lives in a templated C++ routine which is
// exposed to Python through pybind11's NumPy support, avoiding heavier
// dependencies.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <algorithm>
#include <vector>

namespace py = pybind11;

namespace tvd {

// Performs 1-D total variation denoising on `input` using the taut-string
// method described by Condat (2013). The signal is copied by value because the
// algorithm mutates it in-place while constructing the lower and upper
// envelopes. `lambda` tunes the strength of the denoising: larger values yield
// flatter signals.
template <typename T>
std::vector<T> tv1d_denoise(std::vector<T> input, T lambda)
{
    // Number of samples in the input signal
    unsigned int width = input.size();

    // Output buffer and scratch arrays used to track segment starts for the
    // lower and upper envelopes of the taut string.
    std::vector<T> output(width);
    std::vector<unsigned int> indstart_low(width);
    std::vector<unsigned int> indstart_up(width);

    // Indices describing the current state of the envelopes.
    unsigned int j_low = 0, j_up = 0, jseg = 0, indjseg = 0, i = 1, indjseg2, ind;

    // Initial envelope values. `output_low_first` / `output_up_first` represent
    // the first point of the lower/upper envelopes, offset by ±lambda.
    T output_low_first = input[0] - lambda;
    T output_low_curr = output_low_first;
    T output_up_first = input[0] + lambda;
    T output_up_curr = output_up_first;
    T twolambda = static_cast<T>(2) * lambda;

    // Trivial case: single-sample signal.
    if (width == 1)
    {
        output = input;
    }
    else
    {
        indstart_low[0] = 0;
        indstart_up[0] = 0;
        width--;
        // Main loop: process each sample while maintaining envelope stacks.
        for (; i < width; ++i)
        {
            if (input[i] >= output_low_curr)
            {
                if (input[i] <= output_up_curr)
                {
                    // Input lies between envelopes: only update upper part.
                    output_up_curr += (input[i] - output_up_curr) /
                                      (i - indstart_up[j_up] + 1);
                    output[indjseg] = output_up_first;
                    while ((j_up > jseg) &&
                           (output_up_curr <= output[ind = indstart_up[j_up - 1]]))
                        output_up_curr += (output[ind] - output_up_curr) *
                                          static_cast<T>(indstart_up[j_up--] - ind) /
                                          (i - ind + 1);
                    if (j_up == jseg)
                    {
                        // Upper envelope crossed lower envelope: adjust.
                        while ((output_up_curr <= output_low_first) &&
                               (jseg < j_low))
                        {
                            indjseg2 = indstart_low[++jseg];
                            output_up_curr += (output_up_curr - output_low_first) *
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
                    // Current sample above upper envelope: reset upper stack.
                    output_up_curr = output[i] = input[indstart_up[++j_up] = i];
                    output_low_curr += (input[i] - output_low_curr) /
                                       (i - indstart_low[j_low] + 1);
                    output[indjseg] = output_low_first;
                    while ((j_low > jseg) &&
                           (output_low_curr >= output[ind = indstart_low[j_low - 1]]))
                        output_low_curr += (output[ind] - output_low_curr) *
                                           static_cast<T>(indstart_low[j_low--] - ind) /
                                           (i - ind + 1);
                    if (j_low == jseg)
                    {
                        // Lower envelope may cross the upper one: enforce order.
                        while ((output_low_curr >= output_up_first) &&
                               (jseg < j_up))
                        {
                            indjseg2 = indstart_up[++jseg];
                            output_low_curr += (output_low_curr - output_up_first) *
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
                // Symmetric case: sample below lower envelope.
                output_up_curr += ((output_low_curr = output[i] =
                                    input[indstart_low[++j_low] = i]) -
                                   output_up_curr) /
                                  (i - indstart_up[j_up] + 1);
                output[indjseg] = output_up_first;
                while ((j_up > jseg) &&
                       (output_up_curr <= output[ind = indstart_up[j_up - 1]]))
                    output_up_curr += (output[ind] - output_up_curr) *
                                      static_cast<T>(indstart_up[j_up--] - ind) /
                                      (i - ind + 1);
                if (j_up == jseg)
                {
                    while ((output_up_curr <= output_low_first) &&
                           (jseg < j_low))
                    {
                        indjseg2 = indstart_low[++jseg];
                        output_up_curr += (output_up_curr - output_low_first) *
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

        // Finalisation: decide which envelope the last sample belongs to.
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
                output_low_curr += (output[ind] - output_low_curr) *
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
                           (output_up_curr <= output[ind = indstart_up[j_up - 1]]))
                        output_up_curr += (output[ind] - output_up_curr) *
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

    return output;
}

// Convenience wrapper accepting a NumPy array from Python, forwarding it to
// the denoising routine and returning the result as a new array. Templated so
// we can expose overloads for float32 and float64 while sharing the core logic.
template <typename T>
py::array_t<T> TVD(py::array_t<T, py::array::c_style | py::array::forcecast> tinput,
                   double lambda)
{
    auto buf = tinput.request();
    std::vector<T> input(static_cast<T*>(buf.ptr),
                         static_cast<T*>(buf.ptr) + buf.size);
    auto out = tv1d_denoise(std::move(input), static_cast<T>(lambda));
    py::array_t<T> result(out.size());
    std::copy(out.begin(), out.end(), result.mutable_data());
    return result;
}

// Detrend-denoise-retrend helper. First differences of the signal are denoised
// and then cumulatively summed again to reconstruct a smoothed version of the
// original input while roughly preserving the global trend.
template <typename T>
py::array_t<T> D_TVD_R(py::array_t<T, py::array::c_style | py::array::forcecast> tinput,
                        double lambda)
{
    auto buf = tinput.request();
    int size = buf.size;
    const T* data = static_cast<T*>(buf.ptr);
    std::vector<T> detrend(size), retrend(size);

    // Estimate a baseline by averaging the first few samples. This reduces the
    // influence of potential spikes at the beginning of the signal.
    T base_value = data[0];
    T mean = 0;
    size_t sample = static_cast<size_t>(size / 100);
    if (sample < 1)
        sample = 1;
    else if (sample > 10)
        sample = 10;
    for (size_t k = 0; k < sample; ++k)
    {
        mean += data[k];
    }
    mean /= static_cast<T>(sample);
    base_value = mean;

    // Compute first differences (detrending).
    detrend[0] = base_value;
    for (size_t i = 1; i < static_cast<size_t>(size); ++i)
    {
        detrend[i] = data[i] - data[i - 1];
    }

    // Denoise the differences and cumulatively sum to retrend.
    auto denoised = tv1d_denoise(std::move(detrend), static_cast<T>(lambda));

    for (size_t i = 0; i < static_cast<size_t>(size); ++i)
    {
        if (i == 0)
            retrend[i] = base_value;
        else
            retrend[i] = denoised[i] + retrend[i - 1];
    }

    py::array_t<T> result(size);
    std::copy(retrend.begin(), retrend.end(), result.mutable_data());
    return result;
}

} // namespace tvd

PYBIND11_MODULE(TVDCondat2013, m)
{
    m.doc() = R"pbdoc(
        TVDCondat2013 provides 1-D Total Variation Denoising routines
        based on Condat's 2013 algorithm.
    )pbdoc";

    // Bind overloads for both float64 and float32 arrays, offering a uniform
    // Python interface while letting the C++ template pick the right
    // specialisation.
    m.def("TVD", &tvd::TVD<double>, py::arg("signal").noconvert(),
          py::arg("lambda"), "Total variation denoising for float64 arrays");
    m.def("TVD", &tvd::TVD<float>, py::arg("signal").noconvert(),
          py::arg("lambda"), "Total variation denoising for float32 arrays");

    m.def("D_TVD_R", &tvd::D_TVD_R<double>, py::arg("signal").noconvert(),
          py::arg("lambda"), "Detrend-denoise-retrend for float64 arrays");
    m.def("D_TVD_R", &tvd::D_TVD_R<float>, py::arg("signal").noconvert(),
          py::arg("lambda"), "Detrend-denoise-retrend for float32 arrays");
}

