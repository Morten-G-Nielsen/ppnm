#pragma once
#include "core/vector.h"
#include "core/matrix.h"

namespace pp{
  struct EVD{
    vector eigenvalues;
    matrix eigenvectores;

    EVD(matrix A, const double eps);

    const vector& get_values()const{return eigenvalues;}
    const matrix& get_vectores()const{return eigenvectores;}
  };
}
