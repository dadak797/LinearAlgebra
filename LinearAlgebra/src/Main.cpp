#include "pch.h"

#include "Vector.h"
#include <iostream>
#include "BLASInterface.h"


int main()
{
    Vector<DynamicAllocator<double>> matA(9);

    matA = 1.0, -1.0, 0.0
         ,-1.0,  3.0,-2.0
         , 0.0, -2.0, 3.0;

    int ipiv[3];
    //int* ipiv;

    int res = BLASInterface::SYTRF(LAPACK_ROW_MAJOR, 'U', 3, matA.Data(), 3, ipiv);

    std::cout << "return value: " << res << std::endl;
    std::cout << "ipiv: " << ipiv[0] << std::endl;
    std::cout << "ipiv: " << ipiv[1] << std::endl;
    std::cout << "ipiv: " << ipiv[2] << std::endl;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            std::cout << matA[i * 3 + j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}