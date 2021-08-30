#ifndef TRIM_H
#define TRIM_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>



typedef struct contif_t {
    size_t rname;
    size_t rname_len;
    size_t length;
    size_t offset;
    size_t line_bases;
    size_t line_width;
    } Contig;



typedef struct segment_t {
    const char *qname;
    const char *rname;
    const char *cigar;
    const char *mapq;
    const char *pos_ptr;
    char *seq;
    size_t qname_len;
    size_t rname_len;
    size_t cigar_len;
    size_t mapq_len;
    size_t seq_len;
    int32_t pos; // Max size 2^31 - 1 according to sam specifications
    uint16_t flag; // Max size 2^16 - 1 according to sam specifications
    } Segment;


#endif
