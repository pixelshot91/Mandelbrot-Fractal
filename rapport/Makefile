TEXCC := pdflatex
SRC := rapport.tex
PDF := rapport.pdf

DIST_NAME := PRPA-TP1

.PHONY: all
all: $(PDF)

.PHONY: help
help:
	@echo Available rules:
	@echo - all, $(PDF)
	@echo -e '\t(needs pdflatex)'
	@echo - dist: create a tarball containing the sources
	@echo -e '\t(needs GNU tar)'
	@echo - distcheck: verify that the build works
	@echo - mostlyclean: remove all files except the PDF
	@echo - clean

$(PDF): $(SRC)

%.pdf: %.tex
# Some tex files need to get some metadata through a first,
# unfinished compilation, so the compiler is called twice.
	$(TEXCC) $< && $(TEXCC) $<

.PHONY: dist
dist: $(DIST_NAME).tar.gz

.PHONY: distcheck
distcheck:
	make distclean
	make dist
	tar -xvf $(DIST_NAME).tar.gz
	cd $(DIST_NAME) && make

$(DIST_NAME).tar.gz:
	git ls-files							\
		| xargs tar --gzip --transform "s|^|$(DIST_NAME)/|" -cvf $@

.PHONY: distclean
distclean:
	$(RM) -r $(DIST_NAME) $(DIST_NAME).tar.gz

.PHONY: mostlyclean
mostlyclean: distclean
	$(RM) -r auto
	$(RM) $(SRC:.tex=.aux) $(SRC:.tex=.log) $(SRC:.tex=.out) \
		$(SRC:.tex=.toc)

.PHONY: clean
clean: mostlyclean
	$(RM) $(PDF)
	$(RM) $(DIST_NAME).tar.gz
