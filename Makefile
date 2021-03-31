SUBDIRS = pep pgms verify peputil
BIN_FILES  = pep/pep verify/verify peputil/abcps peputil/cpyobs 
BIN_FILES += scripts/pepobs scripts/pepint scripts/pepsol scripts/lnkfort scripts/lnkforts scripts/unlnkfort scripts/vcomp
TARGET_DIR = ~/bin

.PHONY: all new install uninstall test clean cleanall bigtest


all: new
	for dir in $(SUBDIRS); do \
	make -C $$dir -f makefile; \
	done

new:
	git update-index
	make -C pep/ -f makefile $@

install:
	mkdir -p $(TARGET_DIR) && test -w $(TARGET_DIR)
	cd $(TARGET_DIR)
	for binfile in $(BIN_FILES); do \
	ln -sf $(CURDIR)/$$binfile $(TARGET_DIR); \
	done

uninstall:
	rm -f $(TARGET_DIR)/pep $(TARGET_DIR)/verify $(TARGET_DIR)/abcps $(TARGET_DIR)/cpyobs
	rm -f $(TARGET_DIR)/pepobs $(TARGET_DIR)/pepint $(TARGET_DIR)/pepsol $(TARGET_DIR)/*lnkfort*
	rm -f $(TARGET_DIR)/vcomp
bigtest:
	cd bigtest; ./bigtest

clean:
	/bin/rm -f */*.o
	/bin/rm -f *~
	/bin/rm -f bigtest/fort.* bigtest/nbody.v* bigtest/t*.obs*
	/bin/rm -f bigtest/t*.abc bigtest/t*.out bigtest/t*.verout
	/bin/rm -f bigtest/tmi.*les bigtest/tmi.nbody bigtest/tv1.vko
	/bin/rm -f bigtest/tin.ppr bigtest/ttr.sbody bigtest/tmn.[mn]*
	/bin/rm -f bigtest/tfr.*int bigtest/tfr.pvm* bigtest/tfr.tsne
	/bin/rm -f bigtest/tfl.w0

cleanall: clean
	/bin/rm -f pep/pep
	/bin/rm -f peputil/abcps peputil/convw0 peputil/cplnq peputil/cpynbq \
	peputil/cpynrm peputil/cpyobs peputil/ctchar peputil/intrlv \
	peputil/fort.11 peputil/abc.dummy peputil/bigconv.out \
	peputil/prepmnpt peputil/psrcard
	/bin/rm -f verify/verify

