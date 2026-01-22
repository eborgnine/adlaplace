all:
	R -e "roxygen2::roxygenise(\"adlaplace\", load_code=\"source\")"
