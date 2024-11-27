SUBDIRS = distributions

jagsmod_LTLIBRARIES = Gompertz.la

Gompertz_la_SOURCES = Gompertz.cc

Gompertz_la_CPPFLAGS = -I/usr/local/include/JAGS

Gompertz_la_LIBADD = distributions/Gompertzdist.la
if WINDOWS
Gompertz_la_LIBADD += -ljags-$(JAGS_MAJOR) -ljrmath-0
else
Gompertz_la_LIBADD += -ljags -ljrmath
endif

Gompertz_la_LDFLAGS = -module -avoid-version
if WINDOWS
Gompertz_la_LDFLAGS += -no-undefined
endif
