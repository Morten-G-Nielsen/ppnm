#pragma once
#include <functional>
#include <vector>
#include <tuple>
#include "core/vector.h"
#include "core/matrix.h"


namespace pp{
  std::tuple<vector, matrix> lsfit(std::vector<std::function<double(double)>> fs, vector x, vector y, vector dy);
}
