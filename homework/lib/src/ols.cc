#include <tuple>
#include <vector>
#include <functional>
#include "core/vector.h"
#include "core/matrix.h"
#include "qr.h"

namespace pp{
  std::tuple<vector, matrix> lsfit(std::vector<std::function<double(double)>> fs, vector x, vector y, vector dy){
   int n = y.size();
   int m = fs.size(); 
   vector b(n);

   // Constructs the system from the data
   matrix A(n,m); 
   for(int j = 0; j<n; j++){
     double temp = dy[j];
     b[j] = y[j]/temp;
     for(int i = 0; i<m; i++){
       A(j,i) = fs[i](x[j])/temp;
     }
   }
   // Does the least square fit
   QR sys(A);
   matrix R = sys.R;
   vector coeef = sys.solve(b);
   // Finds R_inv by backsubstiution
   matrix R_inv(m,m);
   for(int k = 0; k<m; k++){
     R_inv(k,k) = 1.0/R(k,k);
     for(int i = k-1; i>=0; i--){
       double sum = 0;
       for(int j = i+1; j<k; j++){
         sum += R(i,j)*R_inv(j,k);
       }
       R_inv(i,k) = -sum/R(i,i);
     }
   }
   matrix sigma = R_inv*R_inv.transpose();
   return std::tuple(coeef, sigma);
  }
}
