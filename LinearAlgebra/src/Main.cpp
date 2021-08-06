#include "pch.h"

#include "Vector.h"
#include <iostream>


int main()
{
    Vector<DynamicAllocator<double>> vec(3);

    vec = 1.0, 2.0, 3.0;
    vec += 3.0f;

    Vector<DynamicAllocator<double>> vec2(3);

    vec2 = 5.0, 3.0, 1.0;

    Vector<DynamicAllocator<double>> vec3(3);

    vec3 = vec + vec2;
    
    for (const auto& value: vec3)
    {
        std::cout << value << std::endl;
    }

    return 0;
}