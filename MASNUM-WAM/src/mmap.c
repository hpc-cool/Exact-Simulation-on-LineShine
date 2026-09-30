#include <kupl.h>
//#include <numa.h>
//#include <numaif.h>
#include <sys/mman.h>
int malloc_with_mmap(void **addr, size_t length_bytes, const int pagesz_mb,
                            const int nodeid, const char *buffer_name)
{
    int page_mask = 0;
    size_t length_align = length_bytes;
    void *mapAddress = NULL;

    // printf("Malloc %s on numa-%d, length=%ld\n",
    //        buffer_name, nodeid, length_align);

    if (pagesz_mb > 0) {
        page_mask = 20 + cus_log2( pagesz_mb );
        size_t align_bytes = pagesz_mb * 1024 * 1024;
        length_align = (length_bytes + align_bytes - 1) / align_bytes * align_bytes;
        mapAddress = mmap(NULL, length_align, PROT_READ | PROT_WRITE,
                          MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB | (page_mask << 26),
                          -1, 0);
    } else {
        mapAddress = mmap(NULL, length_align, PROT_READ | PROT_WRITE,
                          MAP_PRIVATE | MAP_ANONYMOUS,
                          -1, 0);
    }

    if (mapAddress == MAP_FAILED) {
        printf("mmap failed for %s on numa-%d, length=%ld. Exit\n",
               buffer_name, nodeid, length_align);
        return -1;
    }

    nodemask_t nodemask;
    struct bitmask bitmask = {NUMA_NUM_NODES, nodemask.n};
    numa_bitmask_clearall(&bitmask);
    numa_bitmask_setbit(&bitmask, nodeid);

    int ret = mbind(mapAddress, length_align, MPOL_BIND, nodemask.n, NUMA_NUM_NODES, 0);
    if (ret != 0) {
        printf( "mbind failed for %s on numa-%d. Exit\n", buffer_name, nodeid );
        return -1;
    }

    #if defined(KPF_LOG_MEM_USAGE)
    printf("Malloc %12ld MB on FIX[%02d] for %s\n",
           length_align/1024/1024, nodeid, buffer_name);
    #endif

    *addr = mapAddress;
    memset( *addr, 1, length_align );
    return 0;
}
