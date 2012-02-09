#ifndef MAB_H
#define MAB_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

#include "cfg.h"

#ifndef MAB_MAX_LINE_LENGTH
# define MAB_MAX_LINE_LENGTH 1024
#endif
#ifndef MAB_MAX_FORMAT_LENGTH
# define MAB_MAX_FORMAT_LENGTH 128
#endif


typedef struct MabStream_ {
    char *filename;     // base filename
    int mode;           // 0 for read, 1 for write

    int count;          // index of currently opened file
    FILE *fp;

    /* Current binary data block */
    size_t n;
    size_t size;
    char endian;
    char fmt[MAB_MAX_FMT_LENGTH];
} MabStream;

MabStream *mab_open(const char *path, const char *mode);

int mab_close(MabStream *ms);

/* Return the stdio file stream associated with a MabStream. */
FILE *mab_get_file(MabStream *ms);

/* Read the header from a MAB file.  If 'n', 'size', 'endian', and 'format' are
 * non-NULL then the number of elements, size of each element, endianness of
 * the machine, and data type are set.  ('format' should be a pointer to a
 * character array of length at least ABN_MAX_FORMAT_LENGTH.)  If 'opts' is
 * non-NULL, it will be filled with any name-value pairs included in the
 * header.  Returns 0 on success, non-zero if an error occurred.  Leaves the
 * file pointer immediately after the header, ready to mab_read() the binary
 * data. */
int mab_read_header(MabStream *stream, size_t* n, size_t* size, char* endian, char* format, Config opts);




/* Write a data header followed by 'n' binary data elements of type given by
 * 'fmt'.  If 'opts' is non-NULL, include the dictionary of name-value pairs in
 * the header.  Returns the number of bytes of binary data written. */
size_t abn_write(FILE* stream, const void* ptr, size_t n, const char* fmt, Config opts);

/* Write a data header for an array of 'n' binary data elements of type given
 * by 'fmt'.  If 'opts' is non-NULL, include the dictionary of name-value pairs
 * in the header.  Leaves the file pointer immediately after the header, ready
 * to fwrite() the binary data. */
void abn_write_header(FILE *stream, size_t n, const char* fmt, Config opts);

/* Read in binary data from a .abn file.  An appropriately sized array is
 * allocated with malloc(), a pointer to which is returned through 'ptr'.  If
 * 'n', 'size', 'endian', and 'fmt' are given then the number of elements, size
 * of each element, endianness of the machine, and data type are set.  ('fmt'
 * should be a pointer to a character array of length at least
 * ABN_MAX_FORMAT_LENGTH.)  If 'opts' is non-NULL, it will be filled with any
 * name-value pairs included in the header.  Returns 0 on success, non-zero if
 * an error occurred. */
int abn_read(FILE *stream, void** ptr, size_t* n, size_t* size, char* endian, char* fmt, Config opts);

/* Read the header from a .abn file.  If 'n', 'size', 'endian', and 'fmt' are
 * given then the number of elements, size of each element, endianness of the
 * machine, and data type are set.  ('fmt' should be a pointer to a character
 * array of length at least ABN_MAX_FORMAT_LENGTH.)  If 'opts' is non-NULL, it
 * will be filled with any name-value pairs included in the header.  Returns 0
 * on success, non-zero if an error occurred.  Leaves the file pointer
 * immediately after the header, ready to fread() the binary data. */
int abn_read_header(FILE *stream, size_t* n, size_t* size, char* endian, char* fmt, Config opts);

/* Utility routine. Return 'L' if machine is little endian, 'B' if big endian. */
char abn_endian();

/* Utility routine. Return the size (in bytes) of the specified data struct, or
 * a negative number if the format string is invalid. */
int abn_calcsize(const char* fmt);

/* Utility routine. */
void abn_byte_swap(char* ptr, size_t n, const char* fmt);

#ifdef __cplusplus
}
#endif

#endif /* MAB_H */
