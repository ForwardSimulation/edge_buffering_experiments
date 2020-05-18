#pragma once

#include <cstdint>
#include <cstdio>
#include <memory>
#include <tskit.h>
#include "tskit_tools.hpp"

struct BirthData
{
    double left, right;
    tsk_id_t child;
    std::int32_t next;
    BirthData(double l, double r, tsk_id_t c);
};

struct EdgeBuffer
{
    std::vector<std::int32_t> first;
    std::vector<BirthData> births;

    EdgeBuffer(std::size_t num_nodes);
};

using edge_buffer_ptr = std::unique_ptr<EdgeBuffer>;

struct ExistingEdges
{
    tsk_id_t parent;
    std::size_t start, stop;
    ExistingEdges(tsk_id_t p, std::size_t start_, std::size_t stop_);
};

struct temp_edges
// Used for calls to tsk_edge_table_set_columns
// Within tskit, we'd just use an edge table.
{
    std::vector<double> left, right;
    std::vector<tsk_id_t> parent, child;
    temp_edges() : left{}, right{}, parent{}, child{}
    {
    }

    void
    clear()
    {
        left.clear();
        right.clear();
        parent.clear();
        child.clear();
    }

    void
    add_edge(double l, double r, tsk_id_t p, tsk_id_t c)
    {
        if (r <= l)
            {
                std::ostringstream o;
                o << "bad left/right " << l << ' ' << r;
                throw std::invalid_argument(o.str());
            }
        left.push_back(l);
        right.push_back(r);
        parent.push_back(p);
        child.push_back(c);
    }

    std::size_t
    size() const
    {
        if (left.size() != right.size() || left.size() != parent.size()
            || left.size() != child.size())
            {
                throw std::runtime_error("invalid size of temporary edges");
            }
        return left.size();
    }
};

std::int32_t
get_buffer_end(const edge_buffer_ptr& new_edges, std::size_t i);

void buffer_new_edge(tsk_id_t parent, double left, double right, tsk_id_t child,
                     edge_buffer_ptr& new_edges);

std::size_t
buffer_new_edge_at(std::size_t loc, double left, double right, tsk_id_t child,
                   edge_buffer_ptr& new_edges);

void stitch_together_edges(const std::vector<tsk_id_t>& alive_at_last_simplification,
                           double max_time, edge_buffer_ptr& new_edges,
                           temp_edges& edge_liftover, table_collection_ptr& tables);
