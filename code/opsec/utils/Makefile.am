CC = $(MPICC)
CXX = $(MPICXX)

# configure script decides which utilities to build
bin_PROGRAMS = $(UTILITIES) legendre
EXTRA_PROGRAMS = dm2abn gmock hpx2img img2hpx legendre makeradial multipole realmock vagc2abn

AM_CPPFLAGS = -I$(top_srcdir)/src -DHAVE_MPI=1
AM_LDFLAGS = -L$(top_builddir)/src -lklpk

legendre_SOURCES = legendre.cpp

dm2abn_SOURCES = dm2abn.cpp

gmock_SOURCES = gmock.cpp grid.cpp
gmock_CXXFLAGS = $(OPENMP_CFLAGS)
gmock_LDADD = $(FFTW2_LIBS)

# Requires CImg.h header and libpng
hpx2img_SOURCES = hpx2img.cpp
hpx2img_LDADD = $(FITSIO_LIBS) $(LIBPNG)

img2hpx_SOURCES = img2hpx.cpp
img2hpx_LDADD = $(FITSIO_LIBS) $(LIBPNG)

makeradial_SOURCES = makeradial.cpp
makeradial_LDADD = $(FITSIO_LIBS)

multipole_SOURCES = multipole.cpp
multipole_CXXFLAGS = $(OPENMP_CFLAGS)

realmock_SOURCES = realmock.cpp
realmock_LDADD = $(FITSIO_LIBS)

vagc2abn_SOURCES = vagc2abn.cpp
