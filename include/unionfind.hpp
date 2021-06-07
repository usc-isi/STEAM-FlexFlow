/*
 * UnionFind.hpp
 *
 * some old code...
 */

#ifndef UNIONFIND_HPP
#define UNIONFIND_HPP

#include <iostream>
#include <vector>
#include <unordered_map>
#include <functional>
#include <iterator>

template<typename Data, class Hash = std::hash<Data> >  // for generality 
class UnionFind {

public:

    /**
     * Default constructor. Trivially create an empty UnionFind structure.
     */
    explicit UnionFind(): size(0), dmap(), uptree(), sz() {}

    /**
     * Create a union find structure by iterating through beg to end. Each unique
     * element of this iteration will be inserted into the union find structure.
     */
    template <typename Iterator>
    explicit UnionFind(const Iterator & b, const Iterator & e) {
        for (auto itr = b; itr != e; itr++) 
            add(*itr);
    }

    /**
     * Destructor. Currently unused.
     */
    ~UnionFind(){}

    /*
     * Add an element to the union find structure, if it's unique.
     * If it's equivalent is already in the structure, ignore it and return false.
     */
    bool add(const Data & d);

    /*
     * convenient wrapper for adding a set of elements with iterator.
     * Return num of element added.
     */
    template <typename Iterator>
    int add(const Iterator & b, const Iterator & e);

    /*
     * Union two disjoint sizes together. This union is union by size.
     * Return true if success (or two data are already in the same union),
     * false if one of the data to be union-ed is not in the set.
     */
    bool join(const Data & d1, const Data & d2);
    
    /*
     * Convenient wrapper for unioning through two iterators.
     * return number of elements in the union.
     */
    template <typename Iterator>
    int joinRange(const Iterator & d1, const Iterator & d2);

    /*
     * assert if two elements are connected. If so, return true. If any element
     * is not in the structure return false.
     */
    bool connected(const Data & d1, const Data & d2);

    /**
     * Returns the size of this union.
     */
    unsigned int getSize() const { return this->size; }

    /**
     * Overloaded << operator, 
     * output three columns of: id, uptree[id], size[id]
     */
    template<typename D, class H>
    friend std::ostream& operator<<(std::ostream& os, const UnionFind<D, H>& u);

private:
    unsigned int size;       // current number of element in the structure
    // map each data onto an integer, for union and find. The integer is
    // equivalently the "id" of a data entry. Data also need to have a ==
    std::unordered_map<Data, unsigned int, Hash> dmap;
    std::vector<unsigned int> uptree;     // vector of representing the uptree
    std::vector<unsigned int> sz;         // vector of size of each tree

    /*
     * "find" an element by trace to it's root in the uptree. Does path compression
     * on it's way of going up. If the element is not in tree, return -1;
     * other wise returns an int that represents the common root.
     */
    int find(const Data & d);
};


/*
 * Add an element to the union find structure, if it's unique.
 * If it's equivalent is already in the structure, ignore it and return false.
 */
template<typename Data, class Hash>
bool UnionFind<Data, Hash>::add(const Data & d) {

    // assert uniqueness
    if (dmap.find(d) != dmap.end()) 
        return false;

    // insert the element if it's not already in.
    dmap[d] = size;
    uptree.push_back(size);
    sz.push_back(1);
    size++;

    return true;
}

/*
 * convenient wrapper for adding a set of elements with iterator.
 * Return num of element added.
 */
template<typename Data, class Hash>
template<typename Iterator>
int UnionFind<Data, Hash>::add(const Iterator & b, const Iterator & e) {
    int added = 0;
    for (auto itr = b; itr != e; itr++) {
        if (add(*itr))
            added++;
    }
    return added;
}

/*
 * Union two disjoint sets together. This union is union by size.
 * Return true if success (or two data are already in the same union),
 * false if one of the data to be union-ed is not in the set.
 */
template<typename Data, class Hash>
bool UnionFind<Data, Hash>::join(const Data & d1, const Data & d2) {
    // integers representing two data; further opt: union their root 
    int i1 = find(d1), i2 = find(d2);
    // if any was not find, return false
    if (i1 == -1 || i2 == -1)
        return false;
    // same element
    if (i1 == i2) {
        return true;
    }
    
    if (sz[i1] >= sz[i2]) {
        uptree[i2] = i1;
        sz[i1] += sz[i2];
    } 
    else {
        uptree[i1] = uptree[i2];
        sz[i2] += sz[i1];
    }
    return true;
}

/*
 * Convenient wrapper for unioning through two iterators.
 * return number of elements in the union.
 */
template<typename Data, class Hash>
template<typename Iterator>
int UnionFind<Data, Hash>::joinRange(const Iterator & b, const Iterator & e) {

    if (b == e)
        return 0;

    int numIn = 0;          // number of element in the union
    auto itr = b;

    while (++itr != e) {
        if(join((*b), (*itr)))
            numIn++;
    }

    return numIn;
}

/*
 * assert if two elements are connected. If so, return true. If any element
 * is not in the structure return false.
 */
template<typename Data, class Hash>
bool UnionFind<Data, Hash>::connected(const Data & d1, const Data & d2) {
    // integers representing two data's root. if same, return true.
    int i1 = find(d1), i2 = find(d2);
    
    if (i1 == -1 || i2 == -1 || i1 != i2)
        return false;

    return true;
}

/*
 * "find" an element by trace to it's root in the uptree. Does path compression
 * on it's way of going up. If the element is not in tree, return -1;
 * other wise returns an int that represents the common root.
 */
template<typename Data, class Hash>
int UnionFind<Data, Hash>::find(const Data & d) {
    // try to grab the mapped int of this element. If can't find, it's invalid.
    int id;

    try {
        id = dmap.at(d);
    }
    catch (const std::out_of_range& oor) {
        return -1;
    }

    // trace to the root and compress the path
    std::vector<int> path;
    while (uptree[id] != id) {
        path.push_back(id);
        id = uptree[id];
    }
    // does the compression
    for (int i = 0; (i + 1) < path.size(); i++) {
        uptree[path[i]] = id;
        sz[path[i+1]] -= sz[path[i]];
    }

    return id;
}

// output three columns of: id, uptree[id], size[id]
template<typename D, class H>
std::ostream& operator<<(std::ostream& os, const UnionFind<D, H>& u) {  
    for (int i = 0; i < u.size; i++) 
        os << i << '\t' << u.uptree[i] << '\t' << u.sz[i] << std::endl;
    return os;
}
#endif