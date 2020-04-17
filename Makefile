package = cartosphere

version = 1.0.0

tarname = $(package)

distdir = $(tarname)-$(version)

export prefix=/usr/local

all clean install uninstall cartosphere:
	$(MAKE) -C src $@

dist: $(distdir).tar.gz

distcheck: $(distdir).tar.gz
	gzip -cd $+ | tar xvf -
	$(MAKE) -C $(distdir) all check
	$(MAKE) -C $(distdir) prefix=\
	$$(PWD)/$(distdir)/_inst install uninstall
	$(MAKE) -C $(distdir) clean
	rm -rf $(distdir)
	@echo "*** Package $(distdir).tar.gz ready for distribution ***"

$(distdir).tar.gz: FORCE $(distdir)
	tar chof - $(distdir) |\
		gzip -9 -c >$(distdir).tar.gz
	rm -rf $(distdir)

$(distdir):
	mkdir -p $(distdir)/src
	cp Makefile $(distdir)
	cp src/Makefile $(distdir)/src
	cp src/main.c $(distdir)/src

FORCE:
	-rm $(distdir).tar.gz &> /dev/null
	-rm -rf $(distdir) &> /dev/null

.PHONY: FORCE all clean dist distcheck

.PHONY: install uninstall
