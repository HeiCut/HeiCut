 /******************************************************************************
 * union_find_sequential.h
 * *
 * Sequential union-find data structure taken and modified from VieCut.
 * *
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@univie.ac.at>
 * Copyright (C) 2017-2018 Alexander Noe <alexander.noe@univie.ac.at>
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/


#ifndef UNION_FIND_SEQUENTIAL_H
#define UNION_FIND_SEQUENTIAL_H

//#pragma once

#include <vector>

#include "union_find_base.h"

// A simple Union-Find datastructure implementation.
// This is sometimes also caled "disjoint sets datastructure.
class UnionFindSequential : public UnionFindBase
{
public:
    explicit UnionFindSequential(unsigned n) : m_parent(n), m_rank(n), m_n(n)
    {
        for (unsigned i = 0; i < m_parent.size(); i++)
        {
            m_parent[i] = i;
            m_rank[i] = 0;
        }
    }
    inline bool Union(unsigned lhs, unsigned rhs)
    {
        int set_lhs = Find(lhs);
        int set_rhs = Find(rhs);
        if (set_lhs != set_rhs)
        {
            if (m_rank[set_lhs] < m_rank[set_rhs])
            {
                m_parent[set_lhs] = set_rhs;
            }
            else
            {
                m_parent[set_rhs] = set_lhs;
                if (m_rank[set_lhs] == m_rank[set_rhs])
                    m_rank[set_lhs]++;
            }
            --m_n;
            return true;
        }
        return false;
    }

    inline unsigned Find(unsigned element)
    {
        if (m_parent[element] != element)
        {
            unsigned retValue = Find(m_parent[element]);
            m_parent[element] = retValue; // path compression
            return retValue;
        }
        return element;
    }

    // Returns:
    //   The total number of sets.
    inline unsigned n()
    {
        return m_n;
    }

private:
    std::vector<unsigned> m_parent;
    std::vector<unsigned> m_rank;

    // Number of elements in UF data structure.
    unsigned m_n;
};

#endif // end of UNION_FIND_SEQUENTIAL_H