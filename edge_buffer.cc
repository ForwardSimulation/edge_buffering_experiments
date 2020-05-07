#include <vector>
#include <stdexcept>
#include "edge_buffer.hpp"

BirthData::BirthData(double l, double r, tsk_id_t c)
    : left{l}, right{r}, child{c}, next{-1}
{
    if (r <= l)
        {
            throw std::invalid_argument("BirthData: right <= left");
        }
}

EdgeBuffer::EdgeBuffer(std::size_t num_nodes) : first(num_nodes, -1), births{}
{
}
