#include <stdexcept>
#include "tskit_tools.hpp"

table_collection_ptr
make_table_collection_ptr(double sequence_length)
{
    table_collection_ptr rv(new tsk_table_collection_t(),
                            [](tsk_table_collection_t* tables) {
                                tsk_table_collection_free(tables);
                                delete tables;
                            });
    int err = tsk_table_collection_init(rv.get(), 0);
    rv->sequence_length = sequence_length;
    if (err != 0)
        {
            throw std::runtime_error("could not initialize tsk_table_collection");
        }
    return rv;
}
