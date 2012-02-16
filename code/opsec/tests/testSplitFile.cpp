/* testSplitFile
 *
 * Test program for SplitFile class.
 * Returns 0 on success, 1 on failure. */

#include <cstdio>
#include <cstdlib>

#include "CuTest.h"
#include "SplitFile.h"


void TestSplitFile(CuTest* tc) {
    /* Generate random data, write it to file, split it, then read it back */
    size_t nbytes = 4096;
    int bytes_per_file = 999;
    char* randout = (char*) malloc(nbytes*sizeof(char));
    char* randin = (char*) malloc(nbytes*sizeof(char));

    srandom(32);
    for(int i = 0; i < nbytes; i++)
        randout[i] = (char) (random() % 256);

    const char* tmpdir = getenv("TMPDIR");
    if(!tmpdir)
        tmpdir = "/tmp";

    char outname[256];
    snprintf(outname, 256, "%s/fullrandom", tmpdir);
    FILE* fout = fopen(outname, "w");
    fwrite(randout, sizeof(char), nbytes, fout);
    fclose(fout);

    char splitname[256];
    snprintf(splitname, 256, "%s/splitrandom", tmpdir);
    char command[512];
    snprintf(command, 512, "split -d -b %d %s %s", bytes_per_file, outname, splitname);
    int status = system(command);

    snprintf(splitname, 256, "%s/splitrandom00", tmpdir);
    SplitFile sf(splitname, "r");
    ssize_t nread = sf.read(randin, nbytes);
    CuAssertIntEquals(tc, (int) nbytes, (int) nread);
    for(int i = 0; i < nbytes; i++)
        CuAssertIntEquals(tc, (int) randout[i], (int) randin[i]);

    snprintf(command, 512, "rm -f %s/fullrandom %s/splitrandom*", tmpdir, tmpdir);
    status = system(command);
}

CuSuite* GetSuite() {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, TestSplitFile);

    return suite;
}

int main(int argc, char **argv) {
    CuString* output = CuStringNew();
    CuSuite* suite = CuSuiteNew();

    CuSuiteAddSuite(suite, GetSuite());

    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s\n", output->buffer);

    return (suite->failCount == 0) ? 0 : 1;
}
