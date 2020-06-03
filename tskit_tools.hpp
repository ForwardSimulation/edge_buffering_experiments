#pragma once

#include <memory>
#include <functional>
#include <tskit/tables.h>

using table_collection_ptr
    = std::unique_ptr<tsk_table_collection_t, std::function<void(tsk_table_collection_t*)>>;

table_collection_ptr make_table_collection_ptr(double sequence_length);
