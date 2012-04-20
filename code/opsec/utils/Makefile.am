CC = $(MPICC)
CXX = $(MPICXX)

# All utility programs
EXTRA_PROGRAMS = dm2abn gmock hpx2img img2hpx legendre makeradial multipole realmock vagc2abn

# Select utility programs to actually build, based on availability of external libraries
bin_PROGRAMS = dm2abn legendre makeradial multipole realmock vagc2abn
if OPSEC_AM_HAVE_FFTW2
    bin_PROGRAMS += gmock
endif
if OPSEC_AM_HAVE_CIMG
if OPSEC_AM_HAVE_LIBPNG
    bin_PROGRAMS += hpx2img img2hpx
endif
endif

AM_CPPFLAGS = -I$(top_srcdir)/src
AM_LDFLAGS = -L$(top_builddir)/src -lklpk

dm2abn_SOURCES = dm2abn.cpp

legendre_SOURCES = legendre.cpp

makeradial_SOURCES = makeradial.cpp
makeradial_LDADD = $(FITSIO_LIBS)

multipole_SOURCES = multipole.cpp
multipole_CXXFLAGS = $(OPENMP_CFLAGS)

realmock_SOURCES = realmock.cpp
realmock_LDADD = $(FITSIO_LIBS)

vagc2abn_SOURCES = vagc2abn.cpp

# Requires FFTW2
gmock_SOURCES = gmock.cpp grid.cpp
gmock_LDADD = $(FFTW2_LIBS)

# Requires CImg.h header and libpng
hpx2img_SOURCES = hpx2img.cpp
hpx2img_LDADD = $(FITSIO_LIBS) $(LIBPNG)

img2hpx_SOURCES = img2hpx.cpp
img2hpx_LDADD = $(FITSIO_LIBS) $(LIBPNG)
