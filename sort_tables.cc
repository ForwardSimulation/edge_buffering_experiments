#include <vector>
#include <stdexcept>
#include <algorithm>
// We need a really good compiler here.
// First, we are checking for C++ >= C++17.
// If that check passes, we use the (absolutely
// incredible) C++17 macro to test for header existence,
// which saves us headaches in our build system.
#if __cplusplus >= 201703L && __has_include(<execution>)
#include <execution>
#endif
#include <tskit.h>

struct _edge
{
    double left, right;
    tsk_id_t parent, child;

    // NOTE: this constuctor must exist or the TBB back-end
    // for the parallel sort won't compile
    _edge() : left{}, right{}, parent{TSK_NULL}, child{TSK_NULL}
    {
    }

    _edge(double l, double r, tsk_id_t p, tsk_id_t c)
        : left{l}, right{r}, parent{p}, child{c}
    {
    }
};

void
sort_tables(tsk_table_collection_t* tables, bool parallel)
// Re-implementation of the copy/sort/copy
// semantics that tskit implements for an edge table.
// If (full) C++17 is available, then we provide
// the option of sorting using the parallel algorithm
// library.
// The parallel method requires a compiler that is not
// on conda.  More seriously, if you do find GCC9 on conda,
// you should NOT use it unless you also recompile ALL C++
// dependencies with it!  You risk runtime crashes otherwise,
// because C++ kinda stinks that way.
{
    // We need some check here to say "If there are edge
    // metadata, throw an exception", or update this to
    // copy the metadata,  too.
    std::vector<_edge> edges;
    edges.reserve(tables->edges.num_rows);
    for (decltype(tables->edges.num_rows) i = 0; i < tables->edges.num_rows; ++i)
        {
            edges.emplace_back(tables->edges.left[i], tables->edges.right[i],
                               tables->edges.parent[i], tables->edges.child[i]);
        }

    // This is our comparison function.  We cannot define an
    // operator< for _edge because we need to bind the node
    // times, so we have to use a functional method.
    // This is a copy of the cmp from fwdpp.  Only difference
    // is the final time comparison (fwdpp table times go forwards).
    const auto cmp = [&tables](const _edge& lhs, const _edge& rhs) {
        auto tl = tables->nodes.time[lhs.parent];
        auto tr = tables->nodes.time[rhs.parent];
        if (tl == tr)
            {
                if (lhs.parent == rhs.parent)
                    {
                        if (lhs.child == rhs.child)
                            {
                                return lhs.left < rhs.left;
                            }
                        return lhs.child < rhs.child;
                    }
                return lhs.parent < rhs.parent;
            }
        return tl < tr;
    };

#if __cplusplus >= 201703L && __has_include(<execution>)
    if (parallel == false)
        {
            std::sort(begin(edges), end(edges), cmp);
        }
    else
        {
            // The support for C++17 parallel algorithms first pops up
            // in GCC9.1 or thereabouts.  No idea about clang.
            // I will assume that this just deletes your hard drive
            // on Windows? ;)
            // The GCC implementation is backed by Intel's excellent
            // TBB library, which becomes a run time dependency
            // and so we need -ltbb.
            // Like OMP, etc., this is "greedy" parallelism, and will
            // auto-determine the size of the thread pool.  I have not
            // experimented if the thread pool size can be controlled
            // using the standard TBB API.  Even if it can, that is not
            // portable, as the implementation of the thread pool
            // is up to the compiler vendor.
            // Unlike OMP, TBB does NOT allow the pool size to be set
            // by an environment variable.
            std::sort(std::execution::par, begin(edges), end(edges), cmp);
        }
#else
    // Default to sequential algorithm.
    std::sort(begin(edges), end(edges), cmp);
#endif

    for (std::size_t i = 0; i < edges.size(); ++i)
        {
            tables->edges.left[i] = edges[i].left;
            tables->edges.right[i] = edges[i].right;
            tables->edges.parent[i] = edges[i].parent;
            tables->edges.child[i] = edges[i].child;
        }
}

