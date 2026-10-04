.DEFAULT_GOAL := check

package := $(shell sed -n 's/^Package: //p' DESCRIPTION)
version := $(shell sed -n 's/^Version: //p' DESCRIPTION)
package_tarball := $(package)_$(version).tar.gz
package_files := DESCRIPTION NAMESPACE README.md \
	$(wildcard R/*.R data/* data-raw/* inst/* inst/*/* man/*.Rd \
	  src/*.cpp src/*.h src/Makevars* tests/*.R tests/testthat/*.R \
	  vignettes/*.Rmd)
check_log := $(package).Rcheck/00check.log
vignette_sources := $(wildcard vignettes/*.Rmd)
vignette_outputs := $(patsubst %.Rmd,%.html,$(vignette_sources))


.PHONY: check
check: $(check_log)

.PHONY: build
build: $(package_tarball)

.PHONY: document
document:
	Rscript -e "library(methods); devtools::document();"

.PHONY: install
install: $(package_tarball)
	R CMD INSTALL $(package_tarball)

.PHONY: preview
preview: $(vignette_outputs)

.PHONY: pkgdown
pkgdown:
	Rscript -e "library(methods); pkgdown::build_site();"

.PHONY: check-revdep
check-revdep: $(package_tarball)
	@mkdir -p revdep
	@rm -rf revdep/*.Rcheck revdep/*.tar.gz
	@cp $(package_tarball) revdep/
	Rscript misc/revdep_check.R

$(package_tarball): $(package_files)
	R CMD build .

$(check_log): $(package_tarball)
	R CMD check --as-cran $(package_tarball)

vignettes/%.html: vignettes/%.Rmd
	Rscript -e "rmarkdown::render('$<')"

.PHONY: tags
tags:
	Rscript -e "utils::rtags(path = 'R', ofile = 'TAGS')"

.PHONY: clean
clean:
	@$(RM) -r *~ */*~ *.Rhistory *.tar.gz src/*.so src/*.o \
	*.Rcheck/ *.Rout .\#* *_cache
