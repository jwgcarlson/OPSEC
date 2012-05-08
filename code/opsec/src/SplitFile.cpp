#include <cassert>
#include <cstdio>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "SplitFile.h"

SplitFile::SplitFile(int suflen, bool numsuf) {
    cur_stream = NULL;
    suffix_length = suflen;
    numeric_suffixes = numsuf;
}

SplitFile::SplitFile(const char* filename, const char* mode, int suflen, bool numsuf) {
    cur_stream = NULL;
    open(filename, mode);
    suffix_length = suflen;
    numeric_suffixes = numsuf;
}

SplitFile::~SplitFile() {
    if(isopen())
        fclose(cur_stream);
}

/* Automatic cast to a stdio stream */
SplitFile::operator FILE*() {
    return cur_stream;
}

bool SplitFile::isopen() {
    return (cur_stream != NULL);
}

const std::string& SplitFile::get_base_filename() const {
    return base_filename;
}

const std::string& SplitFile::get_filename() const {
    return cur_filename;
}

long SplitFile::get_filesize() const {
    return cur_size;
}

void SplitFile::open(const char* filename_, const char* mode_) {
    close();
    base_filename = filename_;
    mode = mode_;
    cur_filename = base_filename;
    cur_index = 0;
    cur_stream = fopen(cur_filename.c_str(), mode.c_str());
    if(cur_stream) {
        struct stat s;
        stat(cur_filename.c_str(), &s);
        cur_size = (long) s.st_size;
    }
}

void SplitFile::reopen() {
    open(base_filename.c_str(), mode.c_str());
}

void SplitFile::close() {
    if(cur_stream != NULL)
        fclose(cur_stream);
    cur_size = 0;
    cur_stream = NULL;
    cur_filename = "";
}

size_t SplitFile::skip(size_t count) {
    if(!isopen())
        return 0;
    if(mode.length() < 1 || mode[0] != 'r') {
        /* Not sure what skip() should do in write mode... */
        return 0;
    }

    /* Keep opening new files until we've skipped the requested number of bytes.
     * Most likely this will happen on the first pass. */
    long skipped = 0;
    while(skipped < count) {
        long nleft = count - skipped;   // number of bytes left to skip
        long pos = ftell(cur_stream);   // current file position
        if(pos + nleft > cur_size) {
            /* Need to open the next file */
            skipped += (cur_size - pos);        // the rest of the current file
            open_next();
            if(!isopen())
                /* Couldn't open the next file, need to return early */
                break;
        }
        else {
            fseek(cur_stream, nleft, SEEK_CUR);
            skipped += nleft;
        }
    }
    return (size_t) skipped;
}

/* Read count bytes of data, storing them in the location given by buf.
 * Return the number of bytes successfully read (this may be less than the
 * requested count). */
size_t SplitFile::read(char* buf, size_t count) {
    assert(mode.length() >= 1 && mode[0] == 'r');
    size_t ntotal = 0;      // total number of bytes read so far
    while(ntotal < count) {
        /* If we don't have a valid current file, we failed */
        if(!isopen() || ferror(cur_stream))
            break;

        /* Try to read data from current file */
        size_t nread = fread(buf, 1, count - ntotal, cur_stream);
        ntotal += (size_t) nread;
        buf += nread;

        /* If we're at the end of the current file and still need more data,
         * try opening the next file */
        if(nread == 0 && feof(cur_stream))
            open_next();

        /* Other exit conditions? */
    }

    return ntotal;
}

size_t SplitFile::write(char* buf, size_t count) {
    assert(mode.length() >= 1 && mode[0] == 'w');
    if(!isopen())
        return 0;
    /* TODO: make this fancier, i.e. open a new output file if the current
     * one will exceed a predefined filesize */
    return fwrite(buf, 1, count, cur_stream);
}

void SplitFile::open_next() {
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

    /* Open new file */
    cur_stream = freopen(cur_filename.c_str(), mode.c_str(), cur_stream);
}
