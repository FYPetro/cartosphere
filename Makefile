PACKAGE = cartosphere
VERSION = 0.1.0

TARNAME = $(PACKAGE)
DISTDIR = $(TARNAME)-$(VERSION)

export prefix=/usr/local

export FFTWINC=-I/usr/local/opt/fftw/include
export FFTWLIB=-L/usr/local/opt/fftw/lib -lfftw3
export S2KITINC=-I/usr/local/include/s2kit10
export S2KITOBJ= 
export EIGENINC=-I/usr/local/include/eigen3

all clean install uninstall cartosphere:
	$(MAKE) -C src $@

dist: $(DISTDIR).tar.gz

distcheck: $(DISTDIR).tar.gz
	gzip -cd $+ | tar xvf -
	$(MAKE) -C $(DISTDIR) all check
	$(MAKE) -C $(DISTDIR) prefix=\
	$$(PWD)/$(DISTDIR)/_inst install uninstall
	$(MAKE) -C $(DISTDIR) clean
	rm -rf $(DISTDIR)
	@echo "*** Package $(DISTDIR).tar.gz ready for distribution ***"

$(DISTDIR).tar.gz: FORCE $(DISTDIR)
	tar chof - $(DISTDIR) |\
		gzip -9 -c >$(DISTDIR).tar.gz
	rm -rf $(DISTDIR)

$(DISTDIR):
	mkdir -p $(DISTDIR)/src
	cp Makefile $(DISTDIR)
	cp src/Makefile $(DISTDIR)/src
	cp src/main.c $(DISTDIR)/src

FORCE:
	-rm $(DISTDIR).tar.gz &> /dev/null
	-rm -rf $(DISTDIR) &> /dev/null

.PHONY: FORCE all clean dist distcheck

.PHONY: install uninstall
