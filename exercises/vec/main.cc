#include <iostream>
#include "vec.h"

int main(){
    vec a {1,0,0};
    vec b {0,1,0};
    vec c = cross(a, b);
    std::cout <<c;
    return 0;
}
