FORT = gfortran
FFLAGS = -fno-automatic -w -frecord-marker=8 -m128bit-long-double
B = ../bigtest
G = ../pgms
P = ../pep
U = .
BIN = .

ABCOBJ = $U/abc.o $U/abcopr.o $U/abcsrt.o $U/again.o $U/autoc.o \
$U/bdypma.o $U/bname.o $U/checka.o $U/daytim.o $U/dft2.o $U/disply.o \
$U/fgen.o $U/fixshp.o $U/fnaxe.o $U/four1.o $U/fourg.o $U/ftinc.o \
$U/getsym.o $U/hrmpma.o $U/hstgrm.o $U/looker.o $U/nampma.o $U/obsio.o \
$U/parmo.o $U/phspma.o $U/plaxis.o $U/plnext.o $U/ploter.o $U/pltdum.o \
$U/prntcn.o $U/random.o $U/rbspma.o $U/rdlib.o $U/select.o $U/setlab.o \
$U/show.o $U/sorter.o $U/sptpma.o $U/srtplt.o $U/titel.o $U/trans.o \
$U/wrkarr.o $U/zerpar.o $U/99suicid.o \
$U/gaussian.o \
$P/blkdata4.o $P/charman.o $P/dtchck.o $P/ebcdi.o $P/mvc.o $P/zfill.o
PSPLTOBJ = $U/blkdapp.o $U/today.o $G/aplot.o
PLTOBJ = $G/symbol.o $G/symbvt.o $U/axis.o $U/line.o $U/nice8.o $U/number.o $U/sqale.o
CPLOBJ = $U/cplntm.o $U/99suicid.o $P/blkdata1.o $P/chncnc.o $P/cross.o $P/dot.o $P/ebcdi.o $P/jnitl.o $P/kp2ki.o $P/prodct.o
#
# Create binary files for BIGTEST
bigtest: abcps convw0 cplnq cpynbq cpynrm cpyobs ctchar intrlv prepmnpt psrcard $B/bigtest.intpvm09 $B/bigtest.intpvm13 $B/bigtest.inttfr $B/bigtest.obsllr $B/bigtest.obspvmdl $B/bigtest.obspvmnc $B/bigtest.obstfl $B/bigtest.obstfr $B/bigtest.obstin $B/bigtest.obstop $B/bigtest.tflw0 $B/bigtest.pprsntin $B/bigtest.tfrtsne $B/bigtest.tmnmoon $B/bigtest.tmnmrot $B/bigtest.tmnntp1 $B/bigtest.tmnntp72 $B/bigtest.ttrsbody $B/ctatint.nbody740 $B/nbody.f311 $B/nbody.f443 $B/nbody.f452
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
abcps: $(ABCOBJ) $(PLTOBJ) $(PSPLTOBJ) $U/qsort.o
	$(FORT) $(FFLAGS) -o $(BIN)/$@ $(ABCOBJ) $(PLTOBJ) $(PSPLTOBJ) $U/qsort.o
#
convw0: convw0.f
	$(FORT) $(FFLAGS) -o $(BIN)/$@ convw0.f
#
cplnq: $(CPLOBJ) $U/dumcpl.o
	$(FORT) $(FFLAGS) -o $(BIN)/$@ $(CPLOBJ) $U/dumcpl.o
#
cpynbq: $U/cpynbdy.o $U/dumcpl.o $U/pr2000.o $P/prodct.o
	$(FORT) $(FFLAGS) -o $(BIN)/$@ $U/cpynbdy.o $U/dumcpl.o $U/pr2000.o $P/prodct.o
#
cpynrm: $U/cpynrm.o $U/99suicid.o $U/dumcpl.o $P/blkdata8.o $P/dtchck.o $P/zfill.o
	$(FORT) $(FFLAGS) -o $(BIN)/$@ $U/cpynrm.o $U/99suicid.o $U/dumcpl.o $P/blkdata8.o $P/dtchck.o $P/zfill.o
#
cpyobs: $U/cpyobs.o $U/fixshp.o $U/obsfrm.o $U/obsio.o $U/zerpar.o $P/dtchck.o $P/zfill.o $U/dumcpl.o $U/99suicid.o
	$(FORT) $(FFLAGS) -o $(BIN)/$@ $U/cpyobs.o $U/fixshp.o $U/obsfrm.o $U/obsio.o $U/zerpar.o $P/dtchck.o $P/zfill.o $U/dumcpl.o $U/99suicid.o
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
psrcard: $U/psrcard.o $U/fixshp.o $U/obsio.o $U/today.o $U/99suicid.o $P/charman.o $P/decodi.o $P/dtchck.o $P/julday.o $P/mdyjul.o $P/mvc.o $P/zfill.o $G/lagrng.o
	$(FORT) $(FFLAGS) -o $(BIN)/$@ $U/psrcard.o $U/fixshp.o $U/obsio.o $U/today.o $U/99suicid.o $P/charman.o $P/decodi.o $P/dtchck.o $P/julday.o $P/mdyjul.o $P/mvc.o $P/zfill.o $G/lagrng.o
#
.f.o:
	$(FORT) $(FFLAGS) -c -o $@ $<
# include dependencies
abc.o: $U/burst.inc $U/bursu.inc $U/checks.inc $U/graphs.inc $U/hrmnsc.inc $U/inodtabc.inc $U/misc.inc $U/namen.inc $U/pltfrm.inc $U/pltlab.inc $U/plttyp.inc $U/skip.inc $U/sortv.inc $U/span.inc
abcopr.o: $U/graphs.inc $U/misc.inc
again.o: $U/graphs.inc $U/inodtabc.inc $U/pltfrm.inc
autoc.o: $U/graphs.inc $U/inodtabc.inc $U/misc.inc $U/pltlab.inc $U/span.inc
axis.o: $U/plttyp.inc
blkdapp.o: $U/plttyp.inc
checka.o: $U/checks.inc $U/namen.inc $U/pepglob.inc $U/t2a.inc $U/t3a.inc
cplntm.o: $P/funcon.inc
cpyobs.o: $U/pepglob.inc $U/t1a.inc $U/t2a.inc $U/t3a.inc $U/t4a.inc
daytim.o: $U/inodtabc.inc $U/plttyp.inc
disply.o: $U/graphs.inc $U/indices.inc $U/inodtabc.inc $U/misc.inc $U/pltlab.inc $U/span.inc
fgen.o: $U/indices.inc $U/inodtabc.inc $U/misc.inc $U/pltlab.inc $U/span.inc
fnaxe.o: $U/graphs.inc $U/pltfrm.inc $U/pltlab.inc $U/plttyp.inc
fourg.o: $U/inodtabc.inc
getsym.o: $U/graphs.inc
hstgrm.o: $U/graphs.inc $U/inodtabc.inc $U/pltlab.inc $U/sortv.inc
line.o: $U/plttyp.inc
looker.o: $U/graphs.inc $U/inodtabc.inc $U/misc.inc $U/pepglob.inc $U/pltlab.inc $U/plttyp.inc $U/sortv.inc $U/span.inc $U/t1a.inc $U/t2a.inc $U/t3a.inc
nampma.o: $P/prmnms.inc $U/pepglob.inc $U/t2a.inc $U/t3a.inc
obsfrm.o: $U/pepglob.inc $U/t1a.inc $U/t2a.inc $U/t3a.inc $U/t4a.inc
obsio.o: $U/pepglob.inc $U/t1a.inc $U/t2a.inc $U/t3a.inc $U/t4a.inc
parmo.o: $U/hrmnsc.inc $U/inodtabc.inc $U/misc.inc $U/span.inc
plaxis.o: $U/burst.inc $U/graphs.inc $U/pltfrm.inc $U/pltlab.inc $U/plttyp.inc $U/sortv.inc
plnext.o: $U/graphs.inc $U/misc.inc $U/pltfrm.inc $U/pltlab.inc $U/plttyp.inc
ploter.o: $U/graphs.inc $U/misc.inc $U/pltfrm.inc $U/sortv.inc
prepmnpt.o: $U/pepglob.inc $U/t1a.inc $U/t2a.inc $U/t3a.inc $U/t4a.inc
prntcn.o: $U/burst.inc $U/bursu.inc $U/checks.inc $U/graphs.inc $U/hrmnsc.inc $U/inodtabc.inc $U/misc.inc $U/namen.inc $U/pltfrm.inc $U/skip.inc $U/sortv.inc $U/span.inc
psrcard.o: $U/pepglob.inc $U/t1a.inc $U/t2a.inc $U/t3a.inc $U/t4a.inc
rdlib.o: $U/checks.inc $U/graphs.inc $U/inodtabc.inc $U/misc.inc $U/namen.inc $U/pepglob.inc $U/skip.inc $U/span.inc $U/t1a.inc $U/t2a.inc $U/t3a.inc $U/t4a.inc
select.o: $U/indices.inc $U/misc.inc $U/span.inc
setlab.o: $U/namen.inc
show.o: $U/graphs.inc $U/indices.inc $U/inodtabc.inc $U/misc.inc $U/pltlab.inc $U/span.inc
sorter.o: $U/graphs.inc
srtplt.o: $U/burst.inc $U/graphs.inc $U/misc.inc $U/pltfrm.inc $U/sortv.inc
titel.o: $U/graphs.inc $U/pltfrm.inc $U/plttyp.inc
trans.o: $U/indices.inc $U/misc.inc $U/span.inc
zerpar.o: $U/pepglob.inc $U/t2a.inc $U/t3a.inc $U/t4a.inc
