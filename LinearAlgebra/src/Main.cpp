#include "pch.h"

#include "Vector.h"
#include <iostream>
#include "BLASInterface.h"


int main()
{
    Vector<DynamicAllocator<double>> matK(4);

    matK = 6.0, -3.0
         ,-3.0,  6.0;

    Vector<DynamicAllocator<double>> matM(4);

    matM = 3.0, 0.0
         , 0.0, 3.0;

    Vector<DynamicAllocator<double>> vecW(2);
 
    int res = BLASInterface::SYGV(LAPACK_ROW_MAJOR, 1, 'V', 'L', 2, matK.Data(), 2, matM.Data(), 2, vecW.Data());

    std::cout << "Return: " << res << std::endl;
       
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            std::cout << matK[i * 2 + j] << " ";
        }
        std::cout << std::endl;
    }

    //for (const auto& value : matA)
      //  std::cout << value << std::endl;

    for (const auto& value : vecW)
        std::cout << value << std::endl;

    return 0;
}