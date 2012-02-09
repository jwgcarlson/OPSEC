/* testcfg
 *
 * Test program for Config and cfg_* routines.  Uses the CuTest unit testing
 * framework.  Returns 0 on success, 1 on failure. */

#include <stdio.h>
#include <string.h>

#include "CuTest.h"
#include "cfg.h"


void TestCfgGetSet(CuTest* tc) {
    Config cfg, copy;

    cfg = cfg_new();
    CuAssertPtrNotNull(tc, cfg);

    cfg_set(cfg, "foo", "bar");
    CuAssertTrue(tc, cfg_has_key(cfg, "foo"));
    CuAssertStrEquals(tc, "bar", cfg_get(cfg, "foo"));

    cfg_set_char  (cfg, "c", -12);
    cfg_set_short (cfg, "s", -1234);
    cfg_set_int   (cfg, "i", -123456);
    cfg_set_long  (cfg, "l", -12345678);
    cfg_set_uchar (cfg, "uc", 12);
    cfg_set_ushort(cfg, "us", 1234);
    cfg_set_uint  (cfg, "ui", 123456);
    cfg_set_ulong (cfg, "ul", 12345678);
    cfg_set_float (cfg, "f",  1e-2);
    cfg_set_double(cfg, "d",  1e-234);
    CuAssertTrue(tc, cfg_get_char  (cfg, "c") == -12);
    CuAssertTrue(tc, cfg_get_short (cfg, "s") == -1234);
    CuAssertTrue(tc, cfg_get_int   (cfg, "i") == -123456);
    CuAssertTrue(tc, cfg_get_long  (cfg, "l") == -12345678);
    CuAssertTrue(tc, cfg_get_uchar (cfg, "uc") == 12);
    CuAssertTrue(tc, cfg_get_ushort(cfg, "us") == 1234);
    CuAssertTrue(tc, cfg_get_uint  (cfg, "ui") == 123456);
    CuAssertTrue(tc, cfg_get_ulong (cfg, "ul") == 12345678);
//    CuAssertTrue(tc, cfg_get_float (cfg, "f")  == 1e-2);
    CuAssertTrue(tc, cfg_get_double(cfg, "d")  == 1e-234);

    cfg_set_format(cfg, "fmt", "%c %d %g", 'a', 12, 32.4);
    CuAssertStrEquals(tc, "a 12 32.4", cfg_get(cfg, "fmt"));

    CuAssertTrue(tc, cfg_has_keys(cfg, "foo,c,s,i,l,uc,us,ui,ul,f,d,fmt", ","));
    copy = cfg_new_copy(cfg);
    CuAssertTrue(tc, cfg_has_keys(copy, "foo,c,s,i,l,uc,us,ui,ul,f,d,fmt", ","));
    cfg_destroy(copy);

//    cfg_write(cfg, stdout);
    cfg_destroy(cfg);
}

void TestCfgGetArray(CuTest* tc) {
    double a[4], b[4], c[4], d[5];
    d[4] = 4;
    Config cfg = cfg_new();
    CuAssertPtrNotNull(tc, cfg);

    cfg_set(cfg, "a", "0 1 2 3");
    cfg_set(cfg, "b", "0, 1, 2, 3");
    cfg_set(cfg, "c", "0, 1, 2, 3,");
    cfg_set(cfg, "d", "0, 1, 2, 3,,");
    cfg_get_array_double(cfg, "a", 4, a);
    cfg_get_array_double(cfg, "b", 4, b);
    cfg_get_array_double(cfg, "c", 4, c);
    cfg_get_array_double(cfg, "d", 5, d);

    CuAssertTrue(tc, a[0] == 0);
    CuAssertTrue(tc, a[1] == 1);
    CuAssertTrue(tc, a[2] == 2);
    CuAssertTrue(tc, a[3] == 3);
    CuAssertTrue(tc, b[0] == 0);
    CuAssertTrue(tc, b[1] == 1);
    CuAssertTrue(tc, b[2] == 2);
    CuAssertTrue(tc, b[3] == 3);
    CuAssertTrue(tc, c[0] == 0);
    CuAssertTrue(tc, c[1] == 1);
    CuAssertTrue(tc, c[2] == 2);
    CuAssertTrue(tc, c[3] == 3);
    CuAssertTrue(tc, d[0] == 0);
    CuAssertTrue(tc, d[1] == 1);
    CuAssertTrue(tc, d[2] == 2);
    CuAssertTrue(tc, d[3] == 3);
    CuAssertTrue(tc, d[4] == 4);

    cfg_destroy(cfg);
}

void TestCfgNewSub(CuTest* tc) {
    Config cfg, subcfg;

    cfg = cfg_new();
    CuAssertPtrNotNull(tc, cfg);

    cfg_set(cfg, "foo.a", "a");
    cfg_set(cfg, "foo.b", "b");
    cfg_set(cfg, "foo.c", "c");
    cfg_set(cfg, "food", "wise");
    cfg_set(cfg, "foo.bill", "shakespeare");

    subcfg = cfg_new_sub(cfg, "foo.", 0);
    CuAssertTrue(tc, cfg_has_key(subcfg, "foo.a"));
    CuAssertTrue(tc, cfg_has_key(subcfg, "foo.b"));
    CuAssertTrue(tc, cfg_has_key(subcfg, "foo.c"));
    CuAssertTrue(tc, cfg_has_key(subcfg, "foo.bill"));
    CuAssertTrue(tc, cfg_get_char(subcfg, "foo.a") == 'a');
    CuAssertTrue(tc, cfg_get_char(subcfg, "foo.b") == 'b');
    CuAssertTrue(tc, cfg_get_char(subcfg, "foo.c") == 'c');
    CuAssertTrue(tc, strcmp(cfg_get(subcfg, "foo.bill"), "shakespeare") == 0);
//    cfg_write(subcfg, stdout);
    cfg_destroy(subcfg);

    subcfg = cfg_new_sub(cfg, "foo.", 1);
    CuAssertTrue(tc, cfg_has_key(subcfg, "a"));
    CuAssertTrue(tc, cfg_has_key(subcfg, "b"));
    CuAssertTrue(tc, cfg_has_key(subcfg, "c"));
    CuAssertTrue(tc, cfg_has_key(subcfg, "bill"));
    CuAssertTrue(tc, cfg_get_char(subcfg, "a") == 'a');
    CuAssertTrue(tc, cfg_get_char(subcfg, "b") == 'b');
    CuAssertTrue(tc, cfg_get_char(subcfg, "c") == 'c');
    CuAssertTrue(tc, strcmp(cfg_get(subcfg, "bill"), "shakespeare") == 0);
//    cfg_write(subcfg, stdout);
    cfg_destroy(subcfg);

    cfg_destroy(cfg);
}

void TestCfgReadWrite(CuTest* tc) {
    Config cfg;
    FILE* fp;

    cfg = cfg_new();
    CuAssertTrue(tc, cfg_read_line(cfg, "foo = bar") == CFG_OKAY);
    CuAssertTrue(tc, cfg_has_key(cfg, "foo"));
    CuAssertStrEquals(tc, "bar", cfg_get(cfg, "foo"));

    fp = tmpfile();
    CuAssertPtrNotNull(tc, fp);
    CuAssertTrue(tc, cfg_write(cfg, fp) == CFG_OKAY);
    cfg_destroy(cfg);

    rewind(fp);
    cfg = cfg_new();
    CuAssertTrue(tc, !cfg_has_key(cfg, "foo"));
    CuAssertTrue(tc, cfg_read(cfg, fp) == CFG_OKAY);
    CuAssertTrue(tc, cfg_has_key(cfg, "foo"));
    CuAssertStrEquals(tc, "bar", cfg_get(cfg, "foo"));

    cfg_destroy(cfg);
}

CuSuite* GetSuite() {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, TestCfgGetSet);
    SUITE_ADD_TEST(suite, TestCfgGetArray);
    SUITE_ADD_TEST(suite, TestCfgNewSub);
    SUITE_ADD_TEST(suite, TestCfgReadWrite);

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
