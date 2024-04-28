/* DisjointSets.cpp */

#include "ParallelDisjointSets.h"
#include <mutex>
#include <ostream>
#include <iostream>


namespace ds
{

DisjointSets::DisjointSets(std::size_t size) : mNodes(size), mNbSets(size)
{   
    for (auto i = std::size_t(0); i < mNodes.size(); ++i)
    {
        auto& node = mNodes[i];
        node.parent.store(i);
        node.rank = 0;
        node.size = 1;
    }
}

std::size_t DisjointSets::find(std::size_t y) const
{

    // Find the root
    size_t t = mNodes[y].parent.load();
    while (mNodes[y].parent.load() != y){
        (mNodes[y].parent).compare_exchange_weak(t, mNodes[t].parent.load());
        y = mNodes[t].parent.load();
    }
   
    return y;
}

void DisjointSets::unite(std::size_t x, std::size_t y)
{   while (true) {
        size_t x1 = find(x);
        size_t y1 = find(y);
        if (x1 == y1) return;
        mNodes[x1].lock.lock();
        mNodes[y1].lock.lock();
        if (mNodes[x1].parent.load() == x && mNodes[y1].parent.load() == y) {
            x = x1;
            y = y1;
            break;
        } 
        else {
            mNodes[x1].lock.unlock();
            mNodes[y1].lock.unlock();
        }
        x = x1;
        y = y1;
    }
    
    if (mNodes[x].parent != mNodes[y].parent)
    {
        if (mNodes[x].rank < mNodes[y].rank)
        {
            mNodes[x].parent = mNodes[y].parent.load();
            mNodes[y].size += mNodes[x].size;
        }
        else
        {   
            size_t x_parent = mNodes[x].parent.load();
            mNodes[y].parent = x_parent;
            mNodes[x].size += mNodes[y].size;
            if (mNodes[x].rank == mNodes[y].rank)
                ++mNodes[x].rank;
        }
        --mNbSets; 
    }
    mNodes[x].lock.unlock();
    mNodes[y].lock.unlock();
}

bool DisjointSets::same(std::size_t x, std::size_t y) const
{
    return find(x) == find(y);
}

std::size_t DisjointSets::getSize() const
{
    return mNodes.size();
}

std::size_t DisjointSets::getSetSize(std::size_t x) const
{
    return mNodes[find(x)].size;
}

std::size_t DisjointSets::getNbSets() const
{
    return mNbSets;
}

} // namespace ds
