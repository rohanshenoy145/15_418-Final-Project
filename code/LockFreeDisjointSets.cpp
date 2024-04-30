/* DisjointSets.cpp */

#include "LockFreeDisjointSets.h"
#include <mutex>
#include <ostream>
#include <iostream>


namespace ds
{

DisjointSets::DisjointSets(std::size_t size) : mNodes(size), mNbSets(size)
{   
    for (auto i = std::size_t(0); i < mNodes.size(); ++i)
    {
        mNodes[i] = new Node;
        mNodes[i].load()->parent.store(i);
        mNodes[i].load()->rank = 0;
        mNodes[i].load()->size = 1;

    }
}

std::size_t DisjointSets::find(std::size_t y) const
{

    // Find the root
    // std::cout<<"find "<< y<<std::endl;
    size_t t = mNodes[y].load()->parent.load();
    while (mNodes[y].load()->parent.load() != y){
        (mNodes[y].load()->parent).compare_exchange_weak(t, mNodes[t].load()->parent.load());
        y = mNodes[t].load()->parent.load();
        // std::cout<<"find"<<std::endl;
    }
   
    return y;
}

bool DisjointSets::update_root(size_t x, size_t old_rank, size_t y, size_t new_rank) 
{
    Node *old = mNodes[x];
    if (old->parent.load() != x || old->rank != old_rank) return false;
    std::atomic<Node*> new_node = new Node;
    new_node.load()->parent.store(y);
    new_node.load()->rank = new_rank;
    
    return mNodes[x].compare_exchange_weak(old, new_node,  std::memory_order_release,
                                           std::memory_order_relaxed);
}

void DisjointSets::unite(std::size_t x, std::size_t y)
{   while (true) {
        size_t x1 = find(x);
        size_t y1 = find(y);
        if (x1 == y1) return;
        // std::cout<<"unite"<<std::endl;
        size_t xr = mNodes[x1].load()->rank;
        size_t yr = mNodes[y1].load()->rank;
        // std::cout<<"unite"<<std::endl;
        if (xr > yr || (xr == yr && x1 > y1)) {
            if (!update_root(y1, yr, x1, yr)) continue;
            if (xr == yr) update_root(x1, xr, x1, xr+1);
             --mNbSets; 
            return;
        }
        else {
            if (!update_root(x1, xr, y1, xr)) continue;
            if (xr == yr) update_root(y1, yr, y1, yr+1);
             --mNbSets; 
            return;
        }
        // std::cout<<"unite"<<std::endl;
    }
    
}

bool DisjointSets::same(std::size_t x, std::size_t y) const
{
    while (true) {
        size_t x1 = find(x);
        size_t y1 = find(y);
        if (x1 == y1) return true;
        if (mNodes[x1].load()->parent.load() == x1) return false;
    }
}

std::size_t DisjointSets::getSize() const
{
    return mNodes.size();
}

std::size_t DisjointSets::getSetSize(std::size_t x) const
{
    return mNodes[find(x)].load()->size;
}

std::size_t DisjointSets::getNbSets() const
{
    return mNbSets;
}

} // namespace ds
