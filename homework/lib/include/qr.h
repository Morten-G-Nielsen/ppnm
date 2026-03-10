#pragma once
#include "core/vector.h"
#include "core/matrix.h"

namespace pp{
struct QR{
  matrix Q;
  matrix R;
  bool is_valid = false;
  QR(const matrix& other);

  vector solve(const vector& b)const;
  double det()const;
  matrix inverse()const;
};
}
