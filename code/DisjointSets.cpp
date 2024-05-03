/* DisjointSets.cpp */

#include "DisjointSets.h"
#include <mutex>

std::mutex coarse_lock;

namespace ds
{

DisjointSets::DisjointSets(std::size_t size) : mNodes(size), mNbSets(size)
{
    for (auto i = std::size_t(0); i < mNodes.size(); ++i)
    {
        auto& node = mNodes[i];
        node.parent = i;
        node.rank = 0;
        node.size = 1;
    }
}

std::size_t DisjointSets::find(std::size_t x) const
{
    // coarse_lock.lock();
    // Find the root
    auto y = x;
    while (mNodes[y].parent != y)
        y = mNodes[y].parent;
    // Path compression
    while (mNodes[x].parent != x)
    {
        auto& node = mNodes[x];
        x = node.parent;
        node.parent = y;
    }
    // coarse_lock.unlock();
    return y;
}

void DisjointSets::unite(std::size_t x, std::size_t y)
{
    auto& rootX = mNodes[find(x)];
    auto& rootY = mNodes[find(y)];
    // coarse_lock.lock();
    if (rootX.parent != rootY.parent)
    {
        if (rootX.rank < rootY.rank)
        {
            rootX.parent = rootY.parent;
            rootY.size += rootX.size;
        }
        else
        {
            rootY.parent = rootX.parent;
            rootX.size += rootY.size;
            if (rootX.rank == rootY.rank)
                ++rootX.rank;
        }
        --mNbSets; 
    }
    // coarse_lock.unlock();
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
