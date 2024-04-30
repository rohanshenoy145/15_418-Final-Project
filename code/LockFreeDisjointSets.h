/* DisjointSets.h */

#pragma once

#include <vector>
#include <mutex>
#include <atomic>
#include <memory>


namespace ds
{

class DisjointSets
{
public:
    DisjointSets(std::size_t size);

    std::size_t find(std::size_t x) const;
    void unite(std::size_t x, std::size_t y);
    bool same(std::size_t x, std::size_t y) const;
    std::size_t getSize() const;
    std::size_t getSetSize(std::size_t x) const;
    std::size_t getNbSets() const;

private:
    struct Node
    {
        std::atomic<size_t> parent;
        std::size_t rank;
        std::size_t size;
    };

    // mutable std::vector<Node> mNodes;
    std::vector<std::atomic<Node*>> mNodes;
    std::size_t mNbSets;
    bool update_root(size_t x, size_t old_rank, size_t y, size_t new_rank);
};

} // namespace ds
