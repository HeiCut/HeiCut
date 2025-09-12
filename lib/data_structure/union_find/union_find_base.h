/******************************************************************************
 * union_find_base.h
 * *
 * Base class for the union-find data structures of VieCut.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#ifndef UNION_FIND_BASE_H
#define UNION_FIND_BASE_H

class UnionFindBase
{
public:
    virtual inline bool Union(unsigned lhs, unsigned rhs) = 0;

    virtual inline unsigned Find(unsigned element) = 0;

    virtual inline unsigned n() = 0;

    virtual ~UnionFindBase() {}
};

#endif // end of UNION_FIND_BASE_H