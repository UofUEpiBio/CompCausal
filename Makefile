help:
	@echo "Available targets:"
	@echo "  docs: Generate documentation using roxygen2."
	@echo "  check: Run R CMD check on the package."
	@echo "  install: Install the package locally."

.PHONY: help docs check install

docs:
	Rscript -e "devtools::document()"

check:
	Rscript -e "devtools::check()"

install:
	Rscript -e "devtools::install()"

