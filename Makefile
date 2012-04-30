install:
	cd ..
	R CMD INSTALL  ./ --clean

test:
	Rscript -e "library('testthat'); require(devtools); test('./');"
