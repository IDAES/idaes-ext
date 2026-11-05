# Make it slightly easier to compile external functions and build docs
#   The external functions all require the ASL, so set the ASL_HEADERS environment
#   variable to point to the ASL build and this should be good to go.
include Makefile.in

.PHONY: all general_helmholtz general_helmholtz_clean cubic cubic_clean functions functions_clean dist dist_clean install clean iapws95 iapws95_clean

all: general_helmholtz cubic functions dist

general_helmholtz:
	$(MAKE) -C ./general_helmholtz
	$(MAKE) -C ./general_helmholtz test

general_helmholtz_clean:
	$(MAKE) -C ./general_helmholtz clean

cubic:
	$(MAKE) -C ./cubic

cubic_clean:
	$(MAKE) -C ./cubic clean

functions:
	$(MAKE) -C ./functions

functions_clean:
	$(MAKE) -C ./functions clean

dist:
	$(RM) -rf dist
	$(MKDIR) dist/
	$(MKDIR) dist/param_data
	$(CP) general_helmholtz/*.so dist/
	$(CP) general_helmholtz/param_data/*.json dist/param_data/
	$(CP) general_helmholtz/param_data/*.nl dist/param_data/
	$(CP) general_helmholtz/param_data/*.py dist/param_data/
	$(CP) cubic/*.so dist/
	$(CP) functions/*.so dist/

install:
	idaes lib-directory --create
	$(CP) dist/* `idaes lib-directory`

dist_clean:
	$(RM) -rf dist

clean: general_helmholtz_clean cubic_clean functions_clean dist_clean
