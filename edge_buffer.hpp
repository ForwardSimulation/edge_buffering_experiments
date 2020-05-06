#pragma once

#include <cstdint>
#include <cstdio>
#include <tskit.h>

struct BirthData
{
    double left, right;
    tsk_id_t child;
    BirthData(double l, double r, tsk_id_t c);
};

struct EdgeBuffer
{
    std::vector<std::int32_t> first, next;
    std::vector<BirthData> births;

    EdgeBuffer(std::size_t num_nodes);
};

