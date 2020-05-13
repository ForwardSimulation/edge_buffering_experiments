#include <vector>
#include <algorithm>
#include <execution> // Requires COMPLETE C++17 support
#include <tskit.h>

struct _edge
{
    double left, right;
    tsk_id_t parent, child;
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
{
    std::vector<_edge> edges;
    edges.reserve(tables->edges.num_rows);
    for (decltype(tables->edges.num_rows) i = 0; i < tables->edges.num_rows; ++i)
        {
            edges.emplace_back(tables->edges.left[i], tables->edges.right[i],
                               tables->edges.parent[i], tables->edges.child[i]);
        }

    if (parallel == false)
        {
            std::sort(begin(edges), end(edges),
                      [&tables](const _edge& lhs, const _edge& rhs) {
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
                      });
        }
    else // Requires COMPLETE C++17 support
        {
           std::sort(std::execution::par, begin(edges), end(edges),
                      [&tables](const _edge& lhs, const _edge& rhs) {
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
                      });
        }

    for (std::size_t i = 0; i < edges.size(); ++i)
        {
            tables->edges.left[i] = edges[i].left;
            tables->edges.right[i] = edges[i].right;
            tables->edges.parent[i] = edges[i].parent;
            tables->edges.child[i] = edges[i].child;
        }
}

