MODULENAME = bude 
MODULE_SOURCES = bude_emc.cpp bude_ff.cpp

BUILT_SOURCES = $(srcdir)/bude_wrap.cxx
SWIG_SOURCES = bude.i 

INCLUDES = -I@top_srcdir@/src/mmlib  -I@top_srcdir@/src/mods/bude
pkgpython_PYTHON = bude.py
pkgpyexec_LTLIBRARIES = _bude.la 
_bude_la_SOURCES = $(MODULE_SOURCES) $(srcdir)/bude_wrap.cxx $(SWIG_SOURCES)
_bude_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/src
_bude_la_LDFLAGS = -module
_bude_la_LIBADD = $(top_srcdir)/src/mmlib/libmmlib.la

noinst_LTLIBRARIES = libforcefields.la
libforcefields_la_SOURCES = bude_ff.cpp bude_ff.h

$(srcdir)/bude_wrap.cxx : $(SWIG_SOURCES)
	$(SWIG)  $(SWIG_PYTHON_OPT) -I$(top_srcdir)/src/  -I$(top_srcdir)/src/mmlib   -o $@ $<


