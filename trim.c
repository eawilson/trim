#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <getopt.h>
#include <ctype.h>
#include <errno.h>
#include <sys/mman.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "zlib.h"

#include "trim.h"

const uint16_t UNMAPPED = 0x4;
const uint16_t MATE_UNMAPPED = 0x8;
const uint16_t REVERSE = 0x10;
const uint16_t MATE_REVERSE = 0x20;
const uint16_t READ1 = 0x40;
const uint16_t READ2 = 0x80;
const uint16_t SECONDARY = 0x100;
const uint16_t FILTERED = 0x200;
const uint16_t SUPPPLEMENTARY = 0x800;
const uint16_t NON_PRIMARY = 0x100 | 0x800; // SECONDARY | SUPPPLEMENTARY;
const uint16_t BOTH_UNMAPPED = 0x4 | 0x8; // UNMAPPED | MATE_UNMAPPED;
const uint16_t READXREVERSEXUNMAPPEDX = 0x40 | 0x80 | 0x10 | 0x20 | 0x4 | 0x8; // READ1 | READ2 | REVERSE | MATE_REVERSE | UNMAPPED | MATE_UNMAPPED;
const uint16_t READX = 0x40 | 0x80; // READ1 | READ2;

const char *CONSUMES_REF = "MDN=X";
const char *CONSUMES_READ = "MIS=X";


bool endswith(const char *text, const char *suffix);
char *parse_segment(char *line, Segment *segment);
const char *cigar_op(const char *cigar, const char **op, int32_t *num);


int main (int argc, char **argv) {
    /*
     * 
     */
    const char *output_filename = "-", *reference_filename = NULL, *input_filename = NULL, *reference = NULL, *error_msg = NULL, *stats_filename = "stats.json";
    const char *cigar = NULL, *cigar_end = NULL, *op = NULL, *second_op = NULL, *final_op = NULL;
    char *index_filename = NULL;
    size_t len = 0;
    FILE *input_sam = NULL, *output_sam = NULL, *reference_index = NULL;//, *stats_file = NULL;
    int reference_fp = 0, trim_len = 0, i = 0;
    off_t reference_len = 0;
    Contig *index = NULL, *current_contig = NULL;
    char *line = NULL, *contig_name_buffer = NULL, *delim = NULL;
    size_t line_max_len = 0, index_len = 0, index_max_len = 0, contig_name_buffer_len = 0, contig_name_buffer_max_len = 0, offset = 0;
    Segment segment = {0};
    bool rc = false;
    int32_t num = 0, pos = 0;
    
    
    
//     , *umi_buffer = NULL, **valid_umis = NULL, *lines[2][4] = {NULL}, *umis[2] = {NULL}, ch = '\0';
//     bool interleaved = false, whitespace = true, validated_umis = true, invalid_short_read = false, invalid_n_read = false;
//     int umi_length = 0, umi_stem_length = 0, min_read_length = 50;
//     int lineno = 0, readno = 0, consecutive_ns = 0, max_consecutive_ns = 2, i = 0;
//     int best = 0, nextbest = 0, edit_distance = 0, j = 0;
//     size_t umi_buffer_len = 0, valid_umis_len = 0, current_input_filename = 0;
//     size_t total_reads = 0, invalid_umi_reads[2] = {0}, invalid_short_reads = 0, invalid_n_reads = 0, length = 0;
//     size_t line_lens[2][4] = {0}, seq_len = 0, qual_len = 0;
//     gzFile output_fastq = NULL, input_fastqs[2] = {NULL};

    // variable needed by strtol
    char *endptr = NULL;
    long val = 0;
    // variables needed by getopt_long
    int option_index = 0, c = 0;
    static struct option long_options[] = {{"output", required_argument, 0, 'o'},
                                           {"reference", required_argument, 0, 'r'},
                                           {"stats", required_argument, 0, 's'},
                                           {"help", no_argument, 0, 'h'},
                                           {0, 0, 0, 0}};

    // Parse optional arguments
    while (c != -1) {
        c = getopt_long(argc, argv, "o:r:s:h", long_options, &option_index);

        switch (c) {
            case 'o':
                if ((!endswith(optarg, ".sam")) && (strcmp(optarg, "-") != 0)) {
                    fprintf(stderr, "%s: output file must be of type sam\n", argv[0]);
                    exit(EXIT_FAILURE);
                    }
                output_filename = optarg;
                break;
                
            case 'r':
                if ((!endswith(optarg, ".fna")) && (!endswith(optarg, ".fa")) && (!endswith(optarg, ".fasta"))) {
                    fprintf(stderr, "%s: reference file must be of type fasta, fna or fa\n", argv[0]);
                    exit(EXIT_FAILURE);
                    }
                reference_filename = optarg;
                break;
                
            case 's':
                if (!endswith(optarg, ".json")) {
                    fprintf(stderr, "%s: stats file must be of type json\n", argv[0]);
                    exit(EXIT_FAILURE);
                    }
                stats_filename = optarg;
                break;
                
            case 'h':
                // help message
                fprintf(stderr, "Program: trim\n");
                fprintf(stderr, "Version: 1.0.0\n");
                fprintf(stderr, "Usage:   trim [options] input_sam\n");
                fprintf(stderr, "Options: -o, --output STR              Name of output sam or '-' for stdout\n");
                fprintf(stderr, "                                       (default '-')\n");
                fprintf(stderr, "         -r, --reference STR           Name of reference genome file ('fasta',\n");
                fprintf(stderr, "                                       'fna' or 'fa' file)\n");
                fprintf(stderr, "         -s, --stats STR               Name of output statistics file \n");
                fprintf(stderr, "                                       (default 'stats.json')\n");
                fprintf(stderr, "         -h, --help                    Display this message\n");
                exit(EXIT_SUCCESS);
                
            case '?':
                // unknown option, getopt_long already printed an error message.
                exit(EXIT_FAILURE);
            }
        }
    
    if (argc - optind == 0) {
        fprintf(stderr, "%s: no input file supplied\n", argv[0]);
        exit(EXIT_FAILURE);
        }
    else if (argc - optind > 1) {
        fprintf(stderr, "%s: greater than one input file supplied\n", argv[0]);
        exit(EXIT_FAILURE);
        }
    input_filename = argv[optind];
    
    if (strcmp(input_filename, "-") == 0) {
        if ((input_sam = fdopen(dup(fileno(stdin)), "rt")) == NULL) {
            fprintf(stderr, "%s: unable to open stdin for reading\n", argv[0]);
            exit(EXIT_FAILURE);
            }
        }
    else {
        if ((input_sam = fopen(input_filename, "rt")) == NULL) {
            fprintf(stderr, "%s: unable to open %s for reading\n", argv[0], input_filename);
            exit(EXIT_FAILURE);
            }
        }
    

    if (reference_filename == NULL) {
        fprintf(stderr, "%s: --reference is a required argument\n", argv[0]);
        exit(EXIT_FAILURE);
        }
    if ((reference_fp = open(reference_filename, O_RDONLY)) == -1) {
        fprintf(stderr, "%s: unable to open %s for reading\n", argv[0], input_filename);
        exit(EXIT_FAILURE);
        }
    if ((reference_len = lseek(reference_fp, 0, SEEK_END)) == -1) {
        fprintf(stderr, "%s: unable to lseek within reference file\n", argv[0]);
        exit(EXIT_FAILURE);
        }
     if ((reference = mmap(NULL, reference_len, PROT_READ, MAP_SHARED, reference_fp, 0)) == MAP_FAILED) {
        fprintf(stderr, "%s: unable to memory map reference file.\n", argv[0]);
        exit(EXIT_FAILURE);
        }
    
    
    len = strlen(reference_filename);
    if ((index_filename = malloc(len + 5)) == NULL) {
        fprintf(stderr, "%s: unable to allocate memory for index filename\n", argv[0]);
        exit(EXIT_FAILURE);
        }
    memcpy(index_filename, reference_filename, len);
    memcpy(index_filename + len, ".fai", 5);
    if ((reference_index = fopen(index_filename, "rt")) == NULL) {
        fprintf(stderr, "%s: unable to open %s for reading\n", argv[0], index_filename);
        exit(EXIT_FAILURE);
        }
    free(index_filename);
    
    
    for (index_len = 0;; ++index_len) {
        if (index_len == 0 || index_len >= index_max_len) {
            index_max_len += 10;
            if ((index = realloc(index, index_max_len * sizeof(Contig))) == NULL) {
                fprintf(stderr, "%s: unable to allocate memory for index struct\n", argv[0]);
                exit(EXIT_FAILURE);
                }
            }
        
        if (getline(&line, &line_max_len, reference_index) == -1) {
            if (feof(reference_index) && index_len) {
                break;
                }
            fprintf(stderr, "%s: error reading index file\n", argv[0]);
            exit(EXIT_FAILURE);
            }
        //fprintf(stderr, "%s\n", line);
        if ((delim = strchr(line, '\t')) == NULL) {
            fprintf(stderr, "%s: malformed index file\n", argv[0]);
            exit(EXIT_FAILURE);
            }
        len = delim - line;
        while (contig_name_buffer_len + len + 1 > contig_name_buffer_max_len) {
            contig_name_buffer_max_len += 100;
            if ((contig_name_buffer = realloc(contig_name_buffer, contig_name_buffer_max_len)) == NULL) {
                fprintf(stderr, "%s: unable to allocate memory for name buffer\n", argv[0]);
                exit(EXIT_FAILURE);
                }
            }
        index[index_len].rname = contig_name_buffer_len;
        index[index_len].rname_len = len;
        memcpy(contig_name_buffer + contig_name_buffer_len, line, len);
        contig_name_buffer[contig_name_buffer_len + len] = '\0';
        contig_name_buffer_len += len + 1;
        for (i = 0; i < 4; ++i) {
            errno = 0;
            val = strtol(delim + 1, &endptr, 10);
            if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == optarg)) {
                fprintf(stderr, "%s: malformed index file\n", argv[0]);
                exit(EXIT_FAILURE);
                }
            delim = endptr;
            *(&index[index_len].length + i) = val;
            }
        }
    current_contig = index;
    
    if (fclose(reference_index) == EOF) {
        fprintf(stderr, "%s: unable to close index file\n", argv[0]);
        exit(EXIT_FAILURE);
        }    
    
    
    if (strcmp(output_filename, "-") == 0) {
        if ((output_sam= fdopen(dup(fileno(stdout)), "wt")) == NULL) {
            fprintf(stderr, "%s: unable to open stdout for writing\n", argv[0]);
            exit(EXIT_FAILURE);
            }
        }
    else {
        if ((output_sam = fopen(output_filename, "wt")) == NULL) {
            fprintf(stderr, "%s: unable to open %s for writing\n", input_filename, argv[0]);
            exit(EXIT_FAILURE);
            }
        }
   
    
    
    while (true) {
        if (getline(&line, &line_max_len, input_sam) == -1) {
            if (feof(input_sam)) {
                break;
                }
            fprintf(stderr, "%s: error reading sam file\n", argv[0]);
            exit(EXIT_FAILURE);
            }
        
        if (line[0] == '@') {
            if (fprintf(output_sam, "%s", line) < 0) {
                fprintf(stderr, "%s: error writing sam file\n", argv[0]);
                exit(EXIT_FAILURE);
                }
            continue;
            }
        
        if ((error_msg = parse_segment(line, &segment)) != NULL) {
            fprintf(stderr, "%s: %s\n", argv[0], error_msg);
            exit(EXIT_FAILURE);
            }
        
        if (*segment.cigar == '*') {
            if (fprintf(output_sam, "%s", line) < 0) {
                fprintf(stderr, "%s: error writing sam file\n", argv[0]);
                exit(EXIT_FAILURE);
                }
            continue;
            }
        
        rc = (bool)(segment.flag & REVERSE);
        cigar = segment.cigar;

        if (!rc) {
            pos = segment.pos - 2; // -1 to get to first reference base and another -1 to finish with the base before the end
            cigar_end = cigar + segment.cigar_len;
            for (; cigar < cigar_end;) {
                final_op = cigar;
                if ((cigar = cigar_op(cigar, &op, &num)) == NULL) {
                    fprintf(stderr, "%s: malformed cigar string  '%.*s'\n", argv[0], (int)segment.cigar_len, segment.cigar);
                    exit(EXIT_FAILURE);
                    }
                if (strchr(CONSUMES_REF, *op) != NULL) {
                    pos += num;
                    }
                }
            }
        else {
            pos = segment.pos;
            if ((second_op = cigar_op(cigar, &op, &num)) == NULL) {
                fprintf(stderr, "%s: malformed cigar string  '%.*s'\n", argv[0], (int)segment.cigar_len, segment.cigar);
                exit(EXIT_FAILURE);
                }
            }
        
        
        if (*op != 'M' || num < 3) {
            if (fprintf(output_sam, "%s", line) < 0) {
                fprintf(stderr, "%s: error writing sam file\n", argv[0]);
                exit(EXIT_FAILURE);
                }
            continue;
            }
        
        if (segment.rname_len == current_contig->rname_len && memcmp(segment.rname, contig_name_buffer + current_contig->rname, segment.rname_len) != 0) {
            for (current_contig = index; current_contig <  index + index_len; ++current_contig) {
                if (segment.rname_len == current_contig->rname_len && memcmp(segment.rname, contig_name_buffer + current_contig->rname, segment.rname_len) == 0) {
                    break;
                    }
                }
            if (current_contig == index + index_len) {
                fprintf(stderr, "%s: %s is not present in reference\n", argv[0], segment.rname);
                exit(EXIT_FAILURE);
                }
            }
        
        if (pos + 1 > current_contig->length) {
            fprintf(stderr, "%s: sam position greater than reference\n", argv[0]);
            exit(EXIT_FAILURE);
            }
        
        
        trim_len = 0;
        if (!rc) {
            offset = current_contig->offset + (((pos - 1) / current_contig->line_bases) * current_contig->line_width) + ((pos - 1) % current_contig->line_bases);
            if (reference[offset] != segment.seq[segment.seq_len - 2]) {
                trim_len = 2;
                }
            else {
                offset = current_contig->offset + ((pos / current_contig->line_bases) * current_contig->line_width) + (pos % current_contig->line_bases);
                if (reference[offset] != segment.seq[segment.seq_len - 1]) {
                    trim_len = 1;
                    }
                }
            }
        else {
            offset = current_contig->offset + ((pos / current_contig->line_bases) * current_contig->line_width) + (pos % current_contig->line_bases);
            if (reference[offset] != segment.seq[1]) {
                trim_len = 2;
                }
            else {
                offset = current_contig->offset + (((pos - 1) / current_contig->line_bases) * current_contig->line_width) + ((pos - 1) % current_contig->line_bases);
                if (reference[offset] != segment.seq[0]) {
                    trim_len = 1;
                    }
                }
            }
        
        if (!trim_len) {
            if (fprintf(output_sam, "%s", line) < 0) {
                fprintf(stderr, "%s: error writing sam file\n", argv[0]);
                exit(EXIT_FAILURE);
                }
            continue;
            }
        
        
        // trim_len > 0 therefore trim read
        if (!rc) {
            fprintf(output_sam, "%.*s%iM%iS%s", (int)(final_op - segment.qname), segment.qname,
                                                (int)num - trim_len, trim_len,
                                                segment.cigar + segment.cigar_len);
            }
        else {
            fprintf(output_sam, "%.*s%u\t%.*s%iS%iM%s", (int)(segment.pos_ptr - segment.qname), segment.qname,
                                                        (unsigned)(pos + trim_len),
                                                        (int)(segment.mapq_len + 1), segment.mapq,
                                                        trim_len, (int)num - trim_len,
                                                        second_op);
            }
        
        }
    
    
    if (munmap((char *)reference, reference_len) == -1) {
        fprintf(stderr, "%s: error unmapping reference file\n", argv[0]);
        exit(EXIT_FAILURE);
        }
    if (close(reference_fp) == -1) {
        fprintf(stderr, "%s: error closing reference file\n", argv[0]);
        exit(EXIT_FAILURE);
        }
        
    if (fclose(input_sam) == EOF) {
        fprintf(stderr, "%s: error closing %s\n", argv[0], input_filename);
        exit(EXIT_FAILURE);
        }    
    if (fclose(output_sam) == EOF) {
        fprintf(stderr, "%s: error closing %s\n", argv[0], output_filename);
        exit(EXIT_FAILURE);
        }    
        
    
    if (stats_filename) {
        }
//     if ((stats_file = fopen(stats_filename, "r+")) == NULL) {
//         if ((stats_file = fopen(stats_filename, "w")) == NULL) {
//             fprintf(stderr, "Error: Unable to open %s for writing\n", stats_filename);
//             exit(EXIT_FAILURE);
//             }
//         fprintf(stats_file, "{\n");
//         }
//     else {
//         fseek(stats_file, 0, SEEK_END);
//         if (ftell(stats_file) == 0) {
//             fprintf(stats_file, "{}\n");
//             }
//         
//         fseek(stats_file, -1, SEEK_CUR);
//         while (true) {
//             ch = fgetc(stats_file);
//             if (ch != '}' && !isspace(ch)) {
//                 if (ch != '{') {
//                     fprintf(stats_file, ",");
//                     }
//                 break;
//                 }
//             if (fseek(stats_file, -2, SEEK_CUR) == -1) {
//                 fprintf(stderr, "Error: Malformed stats file\n");
//                 exit(EXIT_FAILURE);
//                 }
//             }
//         }
// 
//     fprintf(stats_file, "\n    \"total_reads\": %zu,\n", total_reads);
//     fprintf(stats_file, "    \"invalid_short_reads\": %zu,\n", invalid_short_reads);
//     fprintf(stats_file, "    \"invalid_n_reads\": %zu", invalid_n_reads);
//     if (umi_length) {
//         fprintf(stats_file, ",\n    \"invalid_umi_reads_r1\": %zu,\n", invalid_umi_reads[0]);
//         fprintf(stats_file, "    \"invalid_umi_reads_r2\": %zu", invalid_umi_reads[1]);
//         }
//     fprintf(stats_file, "\n}\n");
// 
//     if (fclose(stats_file) == EOF) {
//         fprintf(stderr, "Error: Unable to close stats file\n");
//         exit(EXIT_FAILURE);
//         }
    free(line);    
    free(index);
    free(contig_name_buffer);
    }



bool endswith(const char *text, const char *suffix) {
    int offset = strlen(text) - strlen(suffix);
    return (offset >= 0 && strcmp(text + offset, suffix) == 0);
    }



char *parse_segment(char *read, Segment *segment) {
    int column = 0;
    const char *start = NULL, *endptr = NULL;
    long val = 0;
    
    for (start = read; *read != '\0'; ++read) {
        if (*read == '\t') {
            switch (++column) {
                case 1: // qname
                    segment->qname = start;
                    if ((segment->qname_len = read - start) == 0) {
                        return "missing qname in sam file";
                        }
                    break;
                case 2: // flag
                    errno = 0;
                    val = strtol(start, (char **)&endptr, 10);
                    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == start)) {
                        return "invalid flag in sam file";
                        }
                    segment->flag = (uint16_t)val;
                    break;
                case 3: // rname
                    segment->rname = start;
                    if ((segment->rname_len = read - start) == 0) {
                        return "missing rname in sam file";
                        }
                    break;
                case 4: // pos
                    segment->pos_ptr = start;
                    errno = 0;
                    val = strtol(start, (char **)&endptr, 10);
                    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == start)) {
                        return "invalid pos in sam file";
                        }
                    segment->pos = (int32_t)val;
                    break;
                case 5: // mapq
                    segment->mapq = start;
                    if ((segment->mapq_len = read - start) == 0) {
                        return "missing mapq in sam file";
                        }
                    break;
                case 6: // cigar
                    segment->cigar = start;
                    if ((segment->cigar_len = read - start) == 0) {
                        return "missing cigar in sam file";
                        }
                    break;
                case 10: // seq
                    segment->seq = (char *)start;
                    if ((segment->seq_len = read - start) == 0) {
                        return "missing sequence in sam file";
                        }
                    return NULL;
                }
            start = read + 1;
            }
        }
    return "truncated sam file";
    }



const char *cigar_op(const char *cigar, const char **op, int32_t *num) {
    long val = 0;
    const char *endptr = NULL;

    errno = 0;
    val = strtol(cigar, (char **)&endptr, 10);
    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == cigar)) {
        return NULL;
        }
    
    *num = (int32_t)val;
    *op = endptr;
    return endptr + 1;
    }



