#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "abn.h"

/* TODO: use types of guaranteed size (i.e. int32 instead of int) */

char abn_endian() {
    int x = 1;
    return *((char*)&x) ? 'L' : 'B';
}

int abn_calcsize(const char* fmt) {
    const char* p = fmt;
    int size = 0;               // calculated size of struct
    int num = 1;                // repeat count
    char numbuf[32];            // buffer for reading repeat count
    int numbuflen = 0;
    while(*p != '\0') {
        /* Ignore whitespace */
        while(isspace(*p))
            p++;

        /* Read (optional) repeat count */
        while(isdigit(*p)) {
            numbuf[numbuflen++] = *p;
            p++;
        }
        if(numbuflen == 0)
            num = 1;
        else {
            numbuf[numbuflen] = '\0';
            num = atoi(numbuf);
        }

        switch(*p) {
            case 'x':
            case 'c':
            case 'b':
            case 'B':
                size += num*sizeof(char);
                break;
//            case '?':
//                size += num*sizeof(bool);
//                break;
            case 'h':
            case 'H':
                size += num*sizeof(short);
                break;
            case 'i':
            case 'I':
                size += num*sizeof(int);
                break;
            case 'l':
            case 'L':
                size += num*sizeof(long);
                break;
            case 'q':
            case 'Q':
                size += num*sizeof(long long);
                break;
            case 'f':
                size += num*sizeof(float);
                break;
            case 'd':
                size += num*sizeof(double);
                break;
            default:
                fprintf(stderr, "abn_calcsize: illegal struct format '%s'\n", fmt);
                return -1;
        }

        num = 1;
        numbuflen = 0;
        p++;
    }

    return size;
}

static void bswap(char* data, int size) {
    int i;
    char c;
    if((size % 2) == 1)
        return;
    for(i = 0; i < size/2; i++) {
        c = data[i];
        data[i] = data[size-i-1];
        data[size-i-1] = c;
    }
}

void abn_byte_swap(char* data, size_t n, const char* fmt) {
    int i, j;
    const char* p;
    int fieldlen;
    int fields[1024];
    int numfields;
    int num = 1;                // repeat count
    char numbuf[32];            // buffer for reading repeat count
    int numbuflen = 0;

    int size = abn_calcsize(fmt);
    if(size < 0) {
        fprintf(stderr, "abn_byte_swap: invalid format string, byte-swapping not performed\n");
        return;
    }

    /* Determine the length of each field in the struct */
    numfields = 0;
    p = fmt;
    while(*p != '\0') {
        /* Ignore whitespace */
        while(isspace(*p))
            p++;

        /* Read (optional) repeat count */
        while(isdigit(*p)) {
            numbuf[numbuflen++] = *p;
            p++;
        }
        if(numbuflen == 0)
            num = 1;
        else {
            numbuf[numbuflen] = '\0';
            num = atoi(numbuf);
        }

        switch(*p) {
            case 'x':
            case 'c':
            case 'b':
            case 'B':
                fieldlen = sizeof(char);
                break;
//                case '?':
//                    fieldlen = sizeof(bool);
//                    break;
            case 'h':
            case 'H':
                fieldlen = sizeof(short);
                break;
            case 'i':
            case 'I':
                fieldlen = sizeof(int);
                break;
            case 'l':
            case 'L':
                fieldlen = sizeof(long);
                break;
            case 'q':
            case 'Q':
                fieldlen = sizeof(long long);
                break;
            case 'f':
                fieldlen = sizeof(float);
                break;
            case 'd':
                fieldlen = sizeof(double);
                break;
            default:
                fprintf(stderr, "abn_byte_swap: illegal struct format '%s'\n", fmt);
                return;
        }

        for(j = 0; j < num; j++)
            fields[numfields++] = fieldlen;

        num = 1;
        numbuflen = 0;
        p++;
    }

    /* Iterate through the data array, swapping bytes for each field */
    for(i = 0; i < n; i++) {
        for(j = 0; j < numfields; j++) {
            bswap(data, fields[j]);
            data += fields[j];
        }
    }
}

void abn_write_header(FILE *stream, size_t n, const char* fmt, Config opts) {
    char cleanfmt[ABN_MAX_FORMAT_LENGTH];        // sanitized format string (spaces removed)

    char endian = abn_endian();
    int size;

    const char* src = fmt;
    char* dst = cleanfmt;
    while(*src != '\0') {
        if(!isspace(*src))
            *dst++ = *src;
        src++;
    }
    *dst = '\0';

    size = abn_calcsize(cleanfmt);
    if(size < 0) {
        fprintf(stderr, "abn_write_header: cannot write data header\n");
        return;
    }

    fprintf(stream, "\n");
    if(opts)
        cfg_write(opts, stream);
    fprintf(stream, "#!Begin data block: n = %zd, size = %i, endian = %c, format = %s\n", n, size, endian, cleanfmt);
}

size_t abn_write(FILE* stream, const void* ptr, size_t n, const char* fmt, Config opts) {
    int size = abn_calcsize(fmt);
    size_t nwritten;
    abn_write_header(stream, n, fmt, opts);
    nwritten = fwrite(ptr, size, n, stream);
    return nwritten;
}

int abn_read_header(FILE *stream, size_t* n_, size_t* size_, char* endian_, char* fmt_, Config opts) {
    size_t n, size;
    char endian;
    char fmt[ABN_MAX_FORMAT_LENGTH];
    char line[ABN_MAX_LINE_LENGTH];
    size_t nread = 0;

    while(fgets(line, sizeof(line), stream) != NULL) {
        nread = sscanf(line, "#!Begin data block: n = %zd, size = %zd, endian = %c, format = %s\n", &n, &size, &endian, fmt);
        if(nread == 4)
            break;
        if(opts != NULL)
            cfg_read_line(opts, line);
    }
    
    if(nread != 4) {
        fprintf(stderr, "abn_read_header: no header found\n");
        return -1;
    }

    if(n_) *n_ = n;
    if(size_) *size_ = size;
    if(endian_) *endian_ = endian;
    if(fmt_) strncpy(fmt_, fmt, ABN_MAX_FORMAT_LENGTH);

    return 0;
}

int abn_read(FILE* stream, void** ptr, size_t* n_, size_t* size_, char* endian_, char* fmt_, Config opts) {
    size_t n, size;
    char endian;
    char fmt[ABN_MAX_FORMAT_LENGTH];
    int err;
    size_t nread;

    if(ptr == NULL)
        return 1;

    /* Read the data header */
    err = abn_read_header(stream, &n, &size, &endian, fmt, opts);
    if(err)
        return err;

    /* Allocate memory and read in the data */
    *ptr = malloc(n*size);
    if(*ptr == NULL) {
        fprintf(stderr, "abn_read: failed to allocate %zd bytes of memory\n", n*size);
        perror("system error");
        return 1;
    }
    nread = fread(*ptr, size, n, stream);

    /* Perform byte-swapping if necessary */
    if(endian != abn_endian())
        abn_byte_swap((char*) *ptr, n, fmt);

    if(n_) *n_ = n;
    if(size_) *size_ = size;
    if(endian_) *endian_ = endian;
    if(fmt_) strcpy(fmt_, fmt);

    return 0;
}
