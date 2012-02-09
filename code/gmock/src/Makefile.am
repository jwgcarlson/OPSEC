bin_PROGRAMS = gmock gcorr

gmock_SOURCES =	\
    abn.h array.h   cfg.h grid.h   rng.h sampling.h   Spline.h \
    abn.c array.cpp cfg.c grid.cpp rng.c sampling.cpp Spline.cpp \
    gmock.cpp

gcorr_SOURCES =	\
    abn.h array.h   cfg.h estimator.h   grid.h   rng.h sampling.h   Spline.h \
    abn.c array.cpp cfg.c estimator.cpp grid.cpp rng.c sampling.cpp Spline.cpp \
    gcorr.cpp grid.h particle.h vec3.h


EXTRA_DIST = testrng.c

testrng: testrng.c rng.o
	$(CC) -o testrng testrng.c rng.o -lm
