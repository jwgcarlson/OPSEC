#ifndef ABN_H
#define ABN_H

/* Utility routines for reading and writing Annotated BiNary files (.abn).  An
 * annotated binary file is essentially a raw binary file with a descriptive
 * header in human-readable ASCII.  The goal is that if a competent programmer
 * is given a .abn file with _no_ other documentation, he/she will be able to
 * figure out how to read it into his/her program in a reasonable amount of
 * time.  The file is thus easy to read, while retaining all the advantages of
 * binary data (small file size and machine-precision data).  The header is
 * sufficiently simple and self-explanatory that a raw binary file can be
 * wrapped with a .abn header with almost no thought.
 *
 * The following example shows a valid .abn file, with a brief explanation
 * of the format:

-------------------------------------------------------------------------------
# This is a .abn file.  The header consists of arbitrarily many lines of ASCII
# text.  Blank lines are ignored, and comments start with the '#' character.
# Use comments generously to make your binary files self-documenting!

# You can set variables that help the reading program to interpret the binary
# data.  The following line sets the variable 'foo' to the value 'bar'.
foo = bar
# These key-value pairs are returned as a Config object, which is a dictionary
# of key-value pairs stored as strings.  You can interpret the values as
# various other data types using the cfg_get_*() methods.  More details in
# the cfg.h header.
N = 10
pi = 3.14159

# The header goes on for as long as you want.  The binary data starts
# on the line immediately following the special "#!Begin data block" line
# below.  The parameters n, size, endian, and format on that line indicate
# respectively the number of data elements, the size (in bytes) of each
# element, the endianness in which they are stored ('L' or 'B'), and the
# format of the data.  The format string is interpreted according to the rules
# for the Python 'struct' module.  See the Python documentation for complete
# details.  In the most common case, the binary data consists of an array of
# basic data types, so 'format' will be a single character: 'f' for single
# precision floats, 'd' for doubles, 'i' for ints, etc.

# Okay, enough header, here's the binary data (yes, it's empty!):
#!Begin data block: n = 0, size = 8, endian = L, format = d
-------------------------------------------------------------------------------

 * XXX Finish documentation
 * To read a .abn file into memory...
 * Multiple sections of binary data can be included in a single file...
 * To write a .abn file...
 * Endianness...
 * */

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

#include "cfg.h"

#define ABN_MAX_LINE_LENGTH 1024
#define ABN_MAX_FORMAT_LENGTH 128


/* .abn data struct formats:
 *   See documentation for Python's 'struct' module.
 */

/* Write a data header followed by 'n' binary data elements of type given by
 * 'fmt'.  If 'opts' is non-NULL, include the dictionary of name-value pairs
 * in the header.  Returns the total number of data bytes written. */
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

#endif /* ABN_H */
