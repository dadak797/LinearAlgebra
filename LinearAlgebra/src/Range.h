#pragma once


class Range
{
private:
    int m_First;
    int m_Last;

    enum { TO_END = -1 };

public:
    Range(int first, int last)
        : m_First(first), m_Last(last) {}

    int First() const { return m_First; }
    int Last(int toEnd) const
    {
        if (m_Last == TO_END)
        {
            return toEnd;
        }

        return m_Last;
    }

    int Count(int toEnd) const
    {
        return Last(toEnd) - First() + 1;
    }

    static Range All()
    { 
        return Range(0, TO_END);
    }
};