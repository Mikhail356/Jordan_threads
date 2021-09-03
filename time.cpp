/*Gvozdev 08.11.19*/
#include "jordan.h"

double get_time()
{
    struct rusage buf;
    getrusage(RUSAGE_THREAD,&buf);
    return (double)buf.ru_utime.tv_sec+(double)buf.ru_utime.tv_usec/1000000.;
}

double get_full_time()
{
    struct timeval buf;
    gettimeofday(&buf,0);
    return (double)buf.tv_sec+(double)buf.tv_usec/1000000.;
}
