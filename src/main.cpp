#include "pybind11/pybind11.h"

#include "xtensor/xmath.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xadapt.hpp"

#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pyarray.hpp"
#include "xtensor-python/pyvectorize.hpp"

#include <iostream>
#include <numeric>
#include <cmath>

namespace py = pybind11;
using namespace std;


// Code directly taken from Condat (2013)
xt::pyarray<double> TVD(xt::pyarray<double> tinput,  double lambda) {
  
  unsigned int width = tinput.size();

  vector<double> input(width);
  for (size_t it=0;it<width;it++)
  {
    input[it] = double(tinput[it]);
  }
  if(lambda == 0)
  {
    cout << "Lambda = 0, TVD will only be a minimization by the sum of squares" << endl ;
  }

  vector<double> output(width);
  vector<double> indstart_low(width);
  vector<double> indstart_up(width);
  unsigned int j_low = 0, j_up = 0, jseg = 0, indjseg = 0, i=1, indjseg2, ind;
  double output_low_first = input[0]-lambda;
  double output_low_curr = output_low_first;
  double output_up_first = input[0]+lambda;
  double output_up_curr = output_up_first;
  double twolambda=2.0*lambda;
  if (width==1) 
    {
      output=input;
    }
  else
  {

    indstart_low[0] = 0;
    indstart_up[0] = 0;
    width--;
    for (; i<width; i++) {
        if (input[i]>=output_low_curr) {
          if (input[i]<=output_up_curr) {
              output_up_curr+=(input[i]-output_up_curr)/(i-indstart_up[j_up]+1);
              output[indjseg]=output_up_first;
              while ((j_up>jseg)&&(output_up_curr<=output[ind=indstart_up[j_up-1]]))
                output_up_curr+=(output[ind]-output_up_curr)*
                  ((double)(indstart_up[j_up--]-ind)/(i-ind+1));
              if (j_up==jseg) {
                while ((output_up_curr<=output_low_first)&&(jseg<j_low)) {
                  indjseg2=indstart_low[++jseg];
                output_up_curr+=(output_up_curr-output_low_first)*
                  ((double)(indjseg2-indjseg)/(i-indjseg2+1));
                while (indjseg<indjseg2) output[indjseg++]=output_low_first;
                output_low_first=output[indjseg];
                }
              output_up_first=output_up_curr;
              indstart_up[j_up=jseg]=indjseg;
              } else output[indstart_up[j_up]]=output_up_curr;
          } else
              output_up_curr=output[i]=input[indstart_up[++j_up]=i];
            output_low_curr+=(input[i]-output_low_curr)/(i-indstart_low[j_low]+1);
            output[indjseg]=output_low_first;
            while ((j_low>jseg)&&(output_low_curr>=output[ind=indstart_low[j_low-1]]))
              output_low_curr+=(output[ind]-output_low_curr)*
                  ((double)(indstart_low[j_low--]-ind)/(i-ind+1));
            if (j_low==jseg) {
              while ((output_low_curr>=output_up_first)&&(jseg<j_up)) {
              indjseg2=indstart_up[++jseg];
              output_low_curr+=(output_low_curr-output_up_first)*
                ((double)(indjseg2-indjseg)/(i-indjseg2+1));
              while (indjseg<indjseg2) output[indjseg++]=output_up_first;
              output_up_first=output[indjseg];
              }
              if ((indstart_low[j_low=jseg]=indjseg)==i) output_low_first=output_up_first-twolambda;
              else output_low_first=output_low_curr;
            } else output[indstart_low[j_low]]=output_low_curr;
        } else {
            output_up_curr+=((output_low_curr=output[i]=input[indstart_low[++j_low] = i])-
              output_up_curr)/(i-indstart_up[j_up]+1);
            output[indjseg]=output_up_first;
            while ((j_up>jseg)&&(output_up_curr<=output[ind=indstart_up[j_up-1]]))
              output_up_curr+=(output[ind]-output_up_curr)*
                  ((double)(indstart_up[j_up--]-ind)/(i-ind+1));
            if (j_up==jseg) {
              while ((output_up_curr<=output_low_first)&&(jseg<j_low)) {
              indjseg2=indstart_low[++jseg];
              output_up_curr+=(output_up_curr-output_low_first)*
                ((double)(indjseg2-indjseg)/(i-indjseg2+1));
              while (indjseg<indjseg2) output[indjseg++]=output_low_first;
              output_low_first=output[indjseg];
              }
            if ((indstart_up[j_up=jseg]=indjseg)==i) output_up_first=output_low_first+twolambda;
            else output_up_first=output_up_curr;
            } else output[indstart_up[j_up]]=output_up_curr;
        }
    }
    /* here i==width (with value the actual width minus one) */
    if (input[i]+lambda<=output_low_curr) {
          while (jseg<j_low) {
          indjseg2=indstart_low[++jseg];
          while (indjseg<indjseg2) output[indjseg++]=output_low_first;
          output_low_first=output[indjseg];
      }
      while (indjseg<i) output[indjseg++]=output_low_first;
        output[indjseg]=input[i]+lambda;
    } else if (input[i]-lambda>=output_up_curr) {
      while (jseg<j_up) {
          indjseg2=indstart_up[++jseg];
          while (indjseg<indjseg2) output[indjseg++]=output_up_first;
          output_up_first=output[indjseg];
      }
      while (indjseg<i) output[indjseg++]=output_up_first;
        output[indjseg]=input[i]-lambda;
    } else {
          output_low_curr+=(input[i]+lambda-output_low_curr)/(i-indstart_low[j_low]+1);
          output[indjseg]=output_low_first;
          while ((j_low>jseg)&&(output_low_curr>=output[ind=indstart_low[j_low-1]]))
            output_low_curr+=(output[ind]-output_low_curr)*
                  ((double)(indstart_low[j_low--]-ind)/(i-ind+1));
          if (j_low==jseg) {
            if (output_up_first>=output_low_curr)
              while (indjseg<=i) output[indjseg++]=output_low_curr;
            else {
              output_up_curr+=(input[i]-lambda-output_up_curr)/(i-indstart_up[j_up]+1);
              output[indjseg]=output_up_first;
              while ((j_up>jseg)&&(output_up_curr<=output[ind=indstart_up[j_up-1]]))
                output_up_curr+=(output[ind]-output_up_curr)*
                  ((double)(indstart_up[j_up--]-ind)/(i-ind+1));
              while (jseg<j_up) {
              indjseg2=indstart_up[++jseg];
              while (indjseg<indjseg2) output[indjseg++]=output_up_first;
              output_up_first=output[indjseg];
              }
              indjseg=indstart_up[j_up];
              while (indjseg<=i) output[indjseg++]=output_up_curr;
            }
          } else {
            while (jseg<j_low) {
            indjseg2=indstart_low[++jseg];
            while (indjseg<indjseg2) output[indjseg++]=output_low_first;
            output_low_first=output[indjseg];
            }
            indjseg=indstart_low[j_low];
            while (indjseg<=i) output[indjseg++]=output_low_curr;
          }
    }
  }

  auto output_2 = xt::adapt(output);

  return output_2;
}




// Python Module and Docstrings

PYBIND11_MODULE(TVDCondat2013, m)
{
    xt::import_numpy();

    m.doc() = R"pbdoc(
        TVDCondat2013 is a python portage of the 1D Total Variation Denoising algorithm from Condat 2013: "A Direct Algorithm for 1D Total Variation Denoising" (Sign. Proc. Letters) using xtensor and py11 to bind c++ and numpy.

        .. currentmodule:: TVDCondat2013

        .. autosummary::
           :toctree: _generate

           TVD
    )pbdoc";
    m.def("TVD", TVD);
}
