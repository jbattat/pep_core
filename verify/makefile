FORT = gfortran
FFLAGS = -fno-automatic -w -frecord-marker=8 -m128bit-long-double
#
verify: blkdtav.o cmp.o crd.o dsname.o editl.o iinc.o pcmp.o scmp.o udiff.o vdump.o verf.o verify.o vinput.o ../pep/ebcdi.o
	$(FORT) $(FFLAGS) -o verify blkdtav.o cmp.o crd.o dsname.o editl.o iinc.o pcmp.o scmp.o udiff.o vdump.o verf.o verify.o vinput.o ../pep/ebcdi.o
#
.f.o:
	$(FORT) $(FFLAGS) -c -o $@ $<
#
blkdtav.o: chars.inc dflt.inc ibfl.inc pgfba.inc
cmp.o: lcmp.inc ptrs.inc
crd.o: chars.inc ibfl.inc pgfba.inc ptrs.inc
iinc.o: ibfl.inc ptrs.inc
pcmp.o: lcmp.inc
vdump.o: buffer.inc chars.inc ibfl.inc pgfba.inc ptrs.inc undscr.inc
verf.o: buffer.inc ibfl.inc pgfba.inc ptrs.inc
verify.o: aplst.inc chars.inc ibfl.inc lcmp.inc pgfba.inc ptrs.inc undscr.inc
vinput.o: aplst.inc chars.inc dflt.inc lcmp.inc pgfba.inc ptrs.inc
