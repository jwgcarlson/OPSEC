#ifndef SPLITFILE_H
#define SPLITFILE_H

#include <cstdio>
#include <string>

/* Simple class for reading data from a binary file that may or may not have
 * been broken into pieces, for instance using the Unix utility 'split'. */
class SplitFile {
public:
    SplitFile(int suflen = 2, bool numsuf = true);
    SplitFile(const char* filename, const char* mode, int suflen = 2, bool numsuf = true);
    ~SplitFile();

    /* Automatic cast to a stdio stream */
    operator FILE*();

    const std::string& get_filename() const;
    const std::string& get_base_filename() const;

    bool isopen();
    void open(const char* filename, const char* mode);
    void reopen();
    void close();

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
    int cur_index;
    FILE* cur_stream;
    int suffix_length;
    bool numeric_suffixes;

    /* Try to open the next file in the sequence */
    void open_next();
};

#endif // SPLITFILE_H
