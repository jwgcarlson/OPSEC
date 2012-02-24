#ifndef SPLITFILE_H
#define SPLITFILE_H

#include <cassert>
#include <cstdio>
#include <string>
#include <unistd.h>

/* Simple class for reading data from a binary file that may or may not have
 * been broken into pieces using the Unix utility 'split'. */
class SplitFile {
public:
    SplitFile(int suflen = 2, bool numsuf = true) {
        cur_stream = NULL;
        suffix_length = suflen;
        numeric_suffixes = numsuf;
    }

    SplitFile(const char* filename, const char* mode, int suflen = 2, bool numsuf = true) {
        cur_stream = NULL;
        open(filename, mode);
        suffix_length = suflen;
        numeric_suffixes = numsuf;
    }

    ~SplitFile() {
        if(isopen())
            fclose(cur_stream);
    }

    /* Automatic cast to a stdio stream */
    operator FILE*() {
        return cur_stream;
    }

    bool isopen() {
        return (cur_stream != NULL);
    }

    const std::string& get_filename() const {
        return cur_filename;
    }

    const std::string& get_base_filename() const {
        return base_filename;
    }

    void open(const char* filename_, const char* mode_) {
        close();
        base_filename = filename_;
        mode = mode_;
        cur_filename = base_filename;
        cur_index = 0;
        cur_stream = fopen(cur_filename.c_str(), mode.c_str());
    }

    void reopen() {
        open(base_filename.c_str(), mode.c_str());
    }

    void close() {
        if(cur_stream != NULL)
            fclose(cur_stream);
        cur_stream = NULL;
        cur_filename = "";
    }

    /* Read count bytes of data, storing them in the location given by buf.
     * Return the number of bytes successfully read (this may be less than the
     * requested count). */
    size_t read(char* buf, size_t count) {
        assert(mode.length() >= 1 && mode[0] == 'r');
//        printf("Reading %zd bytes from %s\n", count, base_filename.c_str());
        size_t ntotal = 0;      // total number of bytes read so far
        while(ntotal < count) {
            /* If we don't have a valid current file, we failed */
            if(!isopen() || ferror(cur_stream))
                return ntotal;

            /* Try to read data from current file */
            size_t nread = fread(buf, 1, count - ntotal, cur_stream);
            ntotal += (size_t) nread;
            buf += nread;
//            printf("Read %zd bytes from %s, ntotal = %zd\n", nread, cur_filename.c_str(), ntotal);

            /* If we're at the end of the current file and still need more data,
             * try opening the next file */
            if(nread == 0 && feof(cur_stream))
                open_next();

            /* Other exit conditions? */
        }

        return ntotal;
    }

    /* Write count bytes of data, from the location given by buf.
     * Return the number of bytes successfully written (this may be less than
     * the requested count). */
    size_t write(char* buf, size_t count) {
        assert(mode.length() >= 1 && mode[0] == 'w');
        if(!isopen())
            return 0;
        /* TODO: make this fancier, i.e. open a new output file if the current
         * one will exceed a predefined filesize */
        return fwrite(buf, 1, count, cur_stream);
    }

private:
    std::string base_filename;
    std::string mode;
    std::string cur_filename;
    int cur_index;
    FILE* cur_stream;
    int suffix_length;
    bool numeric_suffixes;

    /* Try to open the next file in the sequence */
    void open_next() {
        if(!isopen())
            return;

        cur_index++;

        /* Generate suffix according to the conventions of 'split'  */
        std::string suffix(suffix_length, '!');
        int base;
        char zero;
        if(numeric_suffixes) {
            base = 10;
            zero = '0';
        }
        else {
            base = 26;
            zero = 'a';
        }
        int x = cur_index;
        for(int i = 0; i < suffix_length; i++) {
            suffix[suffix_length-i-1] = zero + (x % base);
            x /= base;
        }

        /* Substitute suffix into base filename */
        size_t len = cur_filename.length();
        cur_filename.replace(len - suffix_length, suffix_length, suffix);
//        printf("Opening new file, %s\n", cur_filename.c_str());

        /* Open new file */
        cur_stream = freopen(cur_filename.c_str(), mode.c_str(), cur_stream);
    }
};

#endif // SPLITFILE_H
