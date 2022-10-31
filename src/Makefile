# Make it slightly easier to compile external functions and build docs
#   The external functions all require the ASL, so set the ASL_HEADERS environment
#   variable to point to the ASL build and this should be good to go.
include Makefile.in

ALL: general_hemlholtz cubic functions dist

.PHONY: iapws95 iapws95_clean

general_hemlholtz:
	$(MAKE) -C ./general_helmholtz
	$(MAKE) -C ./general_helmholtz test

general_hemlholtz_clean:
	$(MAKE) -C ./general_helmholtz clean

.PHONY: cubic
cubic:
	$(MAKE) -C ./cubic

cubic_clean:
	$(MAKE) -C ./cubic clean

.PHONY: functions
functions:
	$(MAKE) -C ./functions

functions_clean:
	$(MAKE) -C ./functions clean

.PHONY: dist
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

.PHONY: install
install:
	idaes lib-directory --create
	$(CP) dist/* `idaes lib-directory`

dist_clean:
	$(RM) -rf dist

clean: general_hemlholtz_clean cubic_clean functions_clean dist_clean dist_clean
