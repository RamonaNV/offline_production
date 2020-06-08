
# Catch errors from piped commands
SHELL=bash
PF=set -o pipefail;

# When running in multithreaded mode on RHEL 6, GotoBLAS
# occasionally trips a segfault at load time. Try until it works.
$(BASE).single_flux.fits : $(TABLE) fit_flux.py utils.py
	@echo Fitting single-muon flux
	for i in $$(seq 1 100); do \
		$(PF) python fit_flux.py --single $<$(GROUP) $@; \
		exc=$$?; if [[ $$exc -ne 139 ]]; then exit $$exc; fi; done; exit $$exc

$(BASE).bundle_flux.fits : $(TABLE) fit_flux.py utils.py
	@echo Fitting bundle flux
	for i in $$(seq 1 100); do \
		$(PF) python fit_flux.py $<$(GROUP) $@; \
		exc=$$?; if [[ $$exc -ne 139 ]]; then exit $$exc; fi; done; exit $$exc

$(BASE).radius.fits : $(TABLE) fit_rdist.py utils.py
	@echo Fitting bundle radial distribution
	for i in $$(seq 1 100); do \
		$(PF) python fit_rdist.py $<$(GROUP) $@; \
		exc=$$?; if [[ $$exc -ne 139 ]]; then exit $$exc; fi; done; exit $$exc

$(BASE).single_energy.fits : $(TABLE) fit_edist.py utils.py
	@echo Fitting single-muon energy distribution
	for i in $$(seq 1 100); do \
		$(PF) python fit_edist.py --single $<$(GROUP) $@; \
		exc=$$?; if [[ $$exc -ne 139 ]]; then exit $$exc; fi; done; exit $$exc

$(BASE).bundle_energy.fits : $(TABLE) fit_edist.py utils.py
	@echo Fitting bundle energy distribution
	for i in $$(seq 1 100); do \
		$(PF) python fit_edist.py $<$(GROUP) $@; \
		exc=$$?; if [[ $$exc -ne 139 ]]; then exit $$exc; fi; done; exit $$exc

all: $(BASE).single_flux.fits $(BASE).bundle_flux.fits $(BASE).radius.fits $(BASE).single_energy.fits $(BASE).bundle_energy.fits

