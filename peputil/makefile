FORT = gfortran
FFLAGS = -fno-automatic -w -frecord-marker=8 -m128bit-long-double
B = ../bigtest
G = ../pgms
P = ../pep
U = .
BIN = .

CPLOBJ = $U/cplntm.o $U/99suicid.o $U/pr2000.o $U/rotelt.o $P/blkdata1.o $P/chncnc.o $P/cross.o $P/dot.o $P/ebcdi.o $P/jnitl.o $P/kp2ki.o $P/prodct.o
#
# Create binary files for BIGTEST
bigtest: abc.dummy convw0 cplnq cpynbq cpynrm cpyobs ctchar intrlv $B/bigtest.intpvm09 $B/bigtest.intpvm13 $B/bigtest.inttfr $B/bigtest.obsllr $B/bigtest.obspvmdl $B/bigtest.obspvmnc $B/bigtest.obstfl $B/bigtest.obstfr $B/bigtest.obstin $B/bigtest.obstop $B/bigtest.tflw0 $B/bigtest.pprsntin $B/bigtest.tfrtsne $B/bigtest.tmnmoon $B/bigtest.tmnmrot $B/bigtest.tmnntp1 $B/bigtest.tmnntp72 $B/bigtest.ttrsbody $B/ctatint.nbody740 $B/nbody.f311 $B/nbody.f443 $B/nbody.f452
	../bigconv cplnq  $B/bigtest.intpvm09 $B/tfr.pvm09
	../bigconv cplnq  $B/bigtest.intpvm13 $B/tfr.pvm13
	../bigconv cplnq  $B/bigtest.inttfr $B/tfr.int
	../bigconv cpyobs $B/bigtest.obsllr $B/tmn.obsllr
	../bigconv cpyobs $B/bigtest.obspvmdl $B/tfr.obspvmdl
	../bigconv cpyobs $B/bigtest.obspvmnc $B/tfr.obspvmnc
	../bigconv cpyobs $B/bigtest.obstfl $B/tfl.obs
	../bigconv cpyobs $B/bigtest.obstfr $B/tfr.obs
	../bigconv cpyobs $B/bigtest.obstin $B/tin.obs
	../bigconv cpyobs $B/bigtest.obstop $B/top.obstop
	../bigconv convw0 $B/bigtest.tflw0 $B/tfl.w0
	../bigconv cpynrm $B/bigtest.pprsntin $B/tin.ppr
	../bigconv cpynrm $B/bigtest.tfrtsne $B/tfr.tsne
	../bigconv cplnq  $B/bigtest.tmnmoon $B/tmn.moon
	../bigconv cplnq  $B/bigtest.tmnmrot $B/tmn.mrot
	../bigconv cpyobs $B/bigtest.tmnntp1 $B/tmn.ntp1
	../bigconv cpyobs $B/bigtest.tmnntp72 $B/tmn.ntp72
	../bigconv intrlv $B/bigtest.ttrsbody $B/ttr.sbody
	../bigconv ctchar $B/ctatint.nbody740 $B/tfr.ctatint
	../bigconv cpynbq $B/nbody.f311 $B/nbody.v311
	../bigconv cpynbq $B/nbody.f443 $B/nbody.v443
	../bigconv cpynbq $B/nbody.f452 $B/nbody.v452
#
abc.dummy:
	echo echo ABC not implemented > abcps
	chmod +x abcps
	echo> abc.dummy
#
convw0: convw0.f
	$(FORT) $(FFLAGS) -o $(BIN)/$@ convw0.f
#
cplnq: $(CPLOBJ) $U/dumcpl.o
	$(FORT) $(FFLAGS) -o $(BIN)/$@ $(CPLOBJ) $U/dumcpl.o
#
cpynbq: $U/cpynbdy.o $U/dumcpl.o $U/pr2000.o $U/rotelt.o $P/prodct.o
	$(FORT) $(FFLAGS) -o $(BIN)/$@ $U/cpynbdy.o $U/dumcpl.o $U/pr2000.o $U/rotelt.o $P/prodct.o
#
cpynrm: $U/cpynrm.o $U/99suicid.o $U/dumcpl.o $P/blkdata8.o $P/dtchck.o $P/zfill.o
	$(FORT) $(FFLAGS) -o $(BIN)/$@ $U/cpynrm.o $U/99suicid.o $U/dumcpl.o $P/blkdata8.o $P/dtchck.o $P/zfill.o
#
cpyobs: $U/cpyobs.o $U/99suicid.o $U/fixshp.o $U/obsfrm.o $U/obsio.o $U/zerpar.o $P/dtchck.o $P/zfill.o $U/dumcpl.o
	$(FORT) $(FFLAGS) -o $(BIN)/$@ $U/cpyobs.o $U/99suicid.o $U/fixshp.o $U/obsfrm.o $U/obsio.o $U/zerpar.o $P/dtchck.o $P/zfill.o $U/dumcpl.o
#
ctchar: $U/ctchar.o $P/dot.o $P/zfill.o
	$(FORT) $(FFLAGS) -o $(BIN)/$@ $U/ctchar.o $P/dot.o $P/zfill.o
#
intrlv: $U/intrlv.f
	$(FORT) $(FFLAGS) -o $(BIN)/$@ $U/intrlv.f
#
prepmnpt: prepmnpt.f $U/99suicid.o $U/fixshp.o $U/obsio.o $U/today.o $P/dtchck.o $P/julday.o $P/mdyjul.o $P/zfill.o
	$(FORT) $(FFLAGS) -o $(BIN)/$@ prepmnpt.f $U/99suicid.o $U/fixshp.o $U/obsio.o $U/today.o $P/dtchck.o $P/julday.o $P/mdyjul.o $P/zfill.o
#
.f.o:
	$(FORT) $(FFLAGS) -c -o $@ $<
#
cplntm.o: $P/funcon.inc
cpyobs.o: $U/pepglob.inc $U/t1a.inc $U/t2a.inc $U/t3a.inc $U/t4a.inc
obsfrm.o: $U/pepglob.inc $U/t1a.inc $U/t2a.inc $U/t3a.inc $U/t4a.inc
obsio.o: $U/pepglob.inc $U/t1a.inc $U/t2a.inc $U/t3a.inc $U/t4a.inc
prepmnpt.o: $U/pepglob.inc $U/t1a.inc $U/t2a.inc $U/t3a.inc $U/t4a.inc
zerpar.o: $U/pepglob.inc $U/t2a.inc $U/t3a.inc $U/t4a.inc
