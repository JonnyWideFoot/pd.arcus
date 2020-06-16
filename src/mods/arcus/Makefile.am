MODULENAME = arcus 
MODULE_SOURCES = arcus.h 

BUILT_SOURCES = $(srcdir)/arcus_wrap.cxx
SWIG_SOURCES = arcus.i 
SUBDIRS = manipulators workspace filters fileio system forcefields protocols
INCLUDES = -I@top_srcdir@/src/mmlib  -I@top_srcdir@/src/mods/arcus
pkgpython_PYTHON = arcus.py
pkgpyexec_LTLIBRARIES = _arcus.la 
_arcus_la_SOURCES = $(MODULE_SOURCES) $(srcdir)/arcus_wrap.cxx $(SWIG_SOURCES)
_arcus_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/src
_arcus_la_LDFLAGS = -module
_arcus_la_LIBADD = $(top_srcdir)/src/mmlib/libmmlib.la manipulators/libmanipulators.la workspace/libworkspace.la filters/libfilters.la fileio/libfileio.la system/libsystem.la forcefields/libforcefields.la protocols/libprotocols.la

$(srcdir)/arcus_wrap.cxx : $(SWIG_SOURCES)
	$(SWIG)  $(SWIG_PYTHON_OPT) -I$(top_srcdir)/src/  -I$(top_srcdir)/src/mmlib   -o $@ $<
#

#libarcus_la_SOURCES = arcus.h
#libarcus_la_LIBADD =
