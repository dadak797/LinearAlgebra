#pragma once


template<typename AllocatorType>
class Vector;

// For assigning values with comma(,)             
// --- Example ---
// Vector<DynamicAllocator<double>> vector(3);
// vector = 1.0, 2.0, 3.0;

template<typename AllocatorType>
class VectorCommaAssignment
{
private:
    //typedef typename AllocatorType::MemberType DataType;
    using DataType = typename AllocatorType::MemberType;

    Vector<AllocatorType>* m_pVec;
    int m_IP;
    DataType m_Value;

public:
    VectorCommaAssignment(Vector<AllocatorType>* pVec, int ip, const DataType& value)
        : m_pVec(pVec), m_IP(ip), m_Value(value) {}

    ~VectorCommaAssignment()
    {
        if (m_pVec)
        {
            for (int i = m_IP; i < m_pVec->Size(); i++)
            {
                (*m_pVec)(i) = m_Value;
            }
        }
    }

    operator DataType() { return m_Value; }  // Type Casting Operator Overloading

    VectorCommaAssignment operator,(DataType value)
    {
        Vector<AllocatorType>* pVec = m_pVec;
        if (m_IP < m_pVec->Size())
        {
            (*m_pVec)(m_IP) = value;
            m_pVec = nullptr;
        }

        return VectorCommaAssignment(pVec, m_IP + 1, value);
    }
};
