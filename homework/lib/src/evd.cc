#include <cmath>
#include "vector.h"
#include "matrix.h"
#include "evd.h"

namespace pp{
  EVD::EVD(matrix A, const double eps)
    : eigenvalues(A.row_count()), eigenvectores(matrix::identity(A.row_count())){

      if(!A.is_symetric()){
        throw std::runtime_error("Matrix must be symetric for Jacobi EVD");
      }
      int n = A.row_count();
      for(int i = 0; i<A.row_count(); ++i){
        eigenvalues[i] = A(i,i);
      }
      bool changed = true;
      int sweep_count = 0;

      double sum = 0;
      for(int i = 0; i<n; ++i){
        for(int j=i+1; j<n; ++j){
          sum += std::abs(A(i,j));
        }
      }

      while(changed && sweep_count<100){
        changed = false;
        double threshold = 0;
        if(sweep_count<4){
          threshold = (0.2*sum/(n*n))*std::pow(0.5, sweep_count);
        } else if(sweep_count<10){
          threshold = eps*10;
        } else {
          threshold = 0;
        }
        for(int p = 0; p<n-1; ++p){
          double* p_ptr = &eigenvectores(p,0);

          for(int q = p+1; q<n; ++q){
            double apq = A(p,q);
            double abs_apq = std::abs(apq);
            double app = A(p,p);
            double aqq = A(q,q);
            if(abs_apq<threshold) continue;
            if(sweep_count > 4 && (abs_apq <= eps*std::abs(app))
                               && (abs_apq <= eps*std::abs(aqq))){
              A(p,q)=A(q,p)=0.0;
              continue;
            }
            // Performs a rotation if the off-diagnol is larger then tolerance
            if(abs_apq > eps){
              double phi = (aqq-app)/(2*apq);
              double phi_sign = (phi>=0) ? 1.0 : -1.0;
              double t = phi_sign/(std::abs(phi)+std::sqrt(phi*phi+1));
              double c = 1/std::sqrt(t*t+1);
              double s = t*c;
              double tau = s/(1+c);

              // tiny optimization to not repeatly calculate the same values
              double c2 = c*c;
              double s2 = s*s;
              double sc2 = 2*s*c*apq;
              A(p,p) = c2*app + s2*aqq - sc2;
              A(q,q) = s2*app + c2*aqq + sc2;
              A(p,q) = A(q,p) = 0.0;

              // To avoid a if statement that kills performance one loop 
              // is split up into three
              for(int i = 0; i<p; ++i){
                double g = A(i,p);
                double h = A(i,q);
                A(i,p) = g - s*(h+g*tau);
                A(i,q) = h + s*(g-h*tau);
              }
              for(int i = p+1; i<q; ++i){
                double g = A(p,i);
                double h = A(i,q);
                A(p,i) = g - s*(h+g*tau);
                A(i,q) = h + s*(g-h*tau);
              }
              for(int i = q+1; i<n; ++i){
                double g = A(p,i);
                double h = A(q,i);
                A(p,i) = g - s*(h+g*tau);
                A(q,i) = h + s*(g-h*tau);
              }
              // Uses a pointer to read and write to the eigenvectores
              // as a row operation and then transposing it at the end
              double* q_ptr = &eigenvectores(q,0);
              for(int i = 0; i<n; ++i){
                double g = p_ptr[i];
                double h = q_ptr[i];
                p_ptr[i] = g-s*(h+g*tau);
                q_ptr[i] = h+s*(g-h*tau);
              }
              changed = true;
            }
          }
        }
        sweep_count++;
      }
      for(int i = 0; i<n; ++i){
        eigenvalues[i] = A(i,i);
      }
      
      for(int i = 0; i<n-1; i++){
        int k = i;
        double p = eigenvalues[i];
        for(int j = i+1; j<n; j++){
          if(eigenvalues[j]<p){
            k = j;
            p = eigenvalues[j];
          }
        }
        if(k!=i){
          eigenvalues[k] = eigenvalues[i];
          eigenvalues[i] = p;
          for(int col=0; col<n; col++){
            double temp = eigenvectores(i,col);
            eigenvectores(i, col) = eigenvectores(k,col);
            eigenvectores(k,col) = temp;
          }
        }
      }
      // swichece the rows and cols for eigenvectores to get
      // the vectores as columns
      eigenvectores.lazy_transpose();
      //std::cout << "Sweep count = " << sweep_count <<"\n";
    }
}
