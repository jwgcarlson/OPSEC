#ifndef SPLITFILE_H
#define SPLITFILE_H

#include <cstdio>
#include <string>

/* SplitFile
 *
 * Simple class for reading data from a binary file that may or may not have
 * been broken into a series of smaller pieces, for instance using the Unix
 * utility 'split'.
 *
 * The class keeps track of the "current" file within the series, including
 * its file name, file size, its index within the series, and an opened file
 * stream.  When the end of that file is reached, the next file in the series
 * is opened. */
class SplitFile {
public:
    SplitFile(int suflen = 2, bool numsuf = true);
    SplitFile(const char* filename, const char* mode, int suflen = 2, bool numsuf = true);
    ~SplitFile();

    /* Automatic cast to a stdio stream */
    operator FILE*();

    /* Return the name of the first file in the series. */
    const std::string& get_base_filename() const;
    
    /* Return the name of the current file. */
    const std::string& get_filename() const;

    /* Return the size of the current file. */
    long get_filesize() const;

    bool isopen();
    void open(const char* filename, const char* mode);
    void reopen();
    void close();

    /* Skip count bytes of data.  Return the number of bytes actually skipped
     * (this may be less than the requested amount if we reach the end of the
     * file). */
    size_t skip(size_t count);

    /* Read count bytes of data, storing them in the location given by buf.
     * Return the number of bytes successfully read (this may be less than the
     * requested count). */
    size_t read(char* buf, size_t count);

    /* Write count bytes of data, reading them from the location given by buf.
     * Return the number of bytes successfully written (this may be less than
     * the requested count). */
    size_t write(char* buf, size_t count);

private:
    std::string base_filename;
    std::string mode;
    std::string cur_filename;
    long cur_size;
    int cur_index;
    FILE* cur_stream;
    int suffix_length;
    bool numeric_suffixes;

    /* Try to open the next file in the sequence */
    void open_next();
};

#endif // SPLITFILE_H
