#include <vector>
#include "edge_buffer.hpp"

BirthData::BirthData(double l, double r, tsk_id_t c) : left{l}, right{r}, child{c}
{
}

EdgeBuffer::EdgeBuffer(std::size_t num_nodes)
    : first(num_nodes, -1), next(num_nodes, -1), births{}
{
}
