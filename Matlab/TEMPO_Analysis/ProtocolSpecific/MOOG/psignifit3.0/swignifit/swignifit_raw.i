/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */

/* This is the interface file for the swig wrapper to psignifit, swignifit
 */

%module swignifit_raw

%{
#define SWIG_FILE_WITH_INIT
#include "psipp.h"
%}

// custom exception handler
// code need not declare that it will throw an exception
%include "exception.i"
%exception {
    try {
        $action
    // psignift BadArgumentException
    } catch (BadArgumentError& e){
        SWIG_exception(SWIG_ValueError, e.message );
    // All remaining STL exceptions
    } catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    }
}

// make the STL vectors available
%include "std_vector.i"
namespace std {
    %template(vector_double) vector<double>;
    %template(vector_int) vector<int>;
};

// include methods for dealing with double pointers
%include cpointer.i
%pointer_functions(double,doublep);

// This translates BadArgumentError (c++) -> ValueError (python)
// including the error message
// however the function must specify that it throws this error
%typemap(throws) BadArgumentError %{
      PyErr_SetString(PyExc_ValueError, $1.message);
      SWIG_fail;
%}


%include "std_string.i"

// we need to ignore the second constructor for PsiData since swig can't handle
// this type of overloading TODO write a factory method in python that
// implements this functionality
%ignore PsiData::PsiData (std::vector<double> x,
                          std::vector<int>    N,
                          std::vector<double> p,
                          int nAFC);

// We wrap the following headers
%include "data.h"
%include "sigmoid.h"
%include "core.h"
%include "prior.h"
%include "psychometric.h"
%include "optimizer.h"
%include "bootstrap.h"
%include "mcmc.h"
%include "mclist.h"
%include "rng.h"
%include "linalg.h"
%include "getstart.h"
%include "integrate.h"
