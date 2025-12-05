#include<omp.h>

#ifndef FOROMP_HPP
#define FOROMP_HPP
static bool   in_parallel_wrapper()            { return omp_in_parallel() != 0; }
static size_t thread_num_wrapper()             { return static_cast<size_t>(omp_get_thread_num()); }
#endif
