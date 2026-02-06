#include <iostream>
#include <vector>
#include "vec.h"

int main(){
    vec test;
    vec a {1,2,3};
    vec b {4,5,6};

    test.print("test = ");
    a.print("a = ");
    b.print("b = ");
    std::cout << "\n";

    test += a;
    test.print("test += a = ");
    test -= b;
    test.print("test-=b = ");
    test *= 6;
    test /= 3;
    test.print("test = ");
    (-test).print("-test = ");
    test = a - b;
    test.print("test=b-a= ");
    std::cout << "\nnorm(test) = " << norm(test) << "\n";
    std::cout << "dot(test,a) = "<< dot(test,a) << "\n";
    std::cout << "cross(test,a) = "<< cross(test,a) << "\n";
    std::cout << "dot(test, cross(test,a)) = " << dot(test, cross(test,a)) << ", dot(a,cross(test,a)) = "<< dot(a, cross(test,a)) << "\n\n";
    auto basis = GramSchmidt({cross(test,a),a,b});
    std::cout << "basis from GramSchmidt({cross(test,a),a,b})\n";
    for (auto v : basis) std::cout << v << "norm = " << norm(v) <<"\n";
    std::cout << "\ndot(basis[0], basis[1]) = " << dot(basis[0], basis[1]) << "\n";
    std::cout << "dot(basis[0], basis[2]) = " << dot(basis[0], basis[2]) << "\n";
    std::cout << "dot(basis[1], basis[2]) = " << dot(basis[1], basis[2]) << "\n";
    std::cout << "\napprox(a, a*(1+1e-7)) = " << (approx(a, a*(1+1e-7)) ? "true":"false") << "\n";
    std::cout << "approx(a, a*(2)) = " << (approx(a, a*(2)) ? "true":"false") << "\n";
    return 0;
}
