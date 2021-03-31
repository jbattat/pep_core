FORT = gfortran
FFLAGS = -fno-automatic -w -frecord-marker=8
pep: blkdata1.o blkdata2.o blkdata3.o blkdata4.o blkdata7.o blkdata8.o blkdata9.o *.o
	$(FORT) $(FFLAGS) -o pep *.o
	rm -rf newdate
	mv -f blkd1prv blkdata1.f
.f.o:
	$(FORT) $(FFLAGS) -c -o $@ $<
new:
	$(FORT) $(FFLAGS) -c *.f
blkdata1.o: blkdata1.f newdate
newdate:	
	cp blkdata1.f newdate
	rm -rf blkd1prv
	mv blkdata1.f blkd1prv
	sed "s/\(data Lnkdat\).*/\1\/\'`date +%Y%m%d`\'\//" newdate>blkdata1.f

SIMPOBJ = \
a1ut1.o a1ut1f.o a1wwv.o acheck.o adam.o adjast.o adjbdy.o adjfor.o \
adjgrd.o adjhar.o adjphs.o adjrbs.o adjsky.o adjspt.o adjust.o \
alaliz.o analiz.o anaznm.o ancor.o angaz.o angctl.o angdif.o angeq.o \
angmer.o angot.o angvp.o anset.o argmnt.o asbfnc.o asbset.o astcrd.o \
asterp.o astrd1.o atm.o b2rd1.o b2reed.o bcheck.o bdrd1.o bdreed.o \
bdyred.o blkdata1.o blkdata2.o blkdata3.o blkdata7.o bodntl.o bodpun.o \
bodred.o bskip.o calcof.o change.o charman.o check.o chkhms.o chncnc.o \
chncne.o chncnt.o chncor.o cmpar1.o cmpar2.o cmpar3.o cmprst.o \
cmstuff.o cnstrn.o compar.o comrit.o comset.o convrt.o corchn.o \
cpartl.o cross.o ctatf.o ctreed.o ctutf.o d2dms.o decodi.o \
dell.o dguess.o diurn.o diurnp.o dltred.o dms2d.o dnpole.o dot.o \
dotn.o dtchck.o dtpun.o dtrpf.o dumstf16.o eatctl.o eatmdl.o eatmdp.o \
eatmmp.o eatzdl.o eatzmp.o \
ebcdi.o eclprc.o eioctl.o eiondl.o eiondp.o elipt.o emipt.o \
emrd1.o emreed.o emtrp.o endtst.o eqnpun.o eqnred.o eqpcor.o errtot.o \
ertcnt.o ertfn.o ertout.o ertsnt.o eshape.o etutf.o eval.o evtrp.o \
extprc.o flatc.o flats.o formu.o frmbdy.o frmeqn.o frmhar.o frmhed.o \
frmmve.o frmphs.o frmprm.o frmrbs.o frmsit.o frmspt.o frmstr.o \
frstat.o funarg.o funcof.o fwsig.o grdpun.o grdshp.o harntl.o harshp.o \
hdpni2.o hermtf.o horizn.o hrtgcn.o hrtrp.o iau2000a.o imitl.o inital.o initl.o \
inoutv.o input.o inwrit.o itypob.o jbdtst.o jnitl.o julday.o kp2ki.o \
ldeps.o ldpsi.o leg.o legndr.o legnds.o legscl.o lesin.o libchk.o \
libcnd.o linchk.o lmbcor.o lunorb.o lunset.o lvtbdy.o lvthar.o \
lvtprm.o mdyjul.o media.o mnrd1.o mnreed.o mntrp.o mulerr.o mvc.o \
newhed.o newpg.o ninstp.o nmlfix.o nrmadd.o nrmcof.o nrmfrm.o nrmict.o \
nrmrit.o nrmsav.o nrmscl.o nrmset.o numint.o nutat.o nuterp.o nutrm.o \
obsbeg.o obsrd1.o obsred.o optic.o partl.o partl1.o partvl.o patctl.o \
pbdic.o pbdprm.o pcheck.o pcops.o pepmain.o pepnut.o peptic.o peptin.o \
perchk.o perdta.o pexlic.o phspun.o phsred.o planet.o plcmc.o plcms.o plink.o \
plnhed.o plrd1.o plreed.o pltrp.o pnitl.o prces9.o prcnt9.o prcnut.o \
prdbdy.o prdeqs.o prdgrd.o prdhar.o prdict.o prdobs.o prdprm.o \
prdset.o prdsid.o prdvct.o preces.o prmred.o prodct.o propco.o \
prpfrq.o prplog.o prpstr.o prpzen.o prpzro.o prtcrd.o prterb.o \
prterp.o prtfn.o prtout.o prtpr1.o prtrd1.o prtrel.o psrpun.o punch.o \
qread.o radar.o radctl.o radns.o radpl.o radrel.o radvaf.o rbspun.o \
rbsred.o redcpy.o reflcc.o refrct.o restrt.o rngmod.o rotat.o rroad.o \
rtrd1.o rtreed.o savhed.o sbeng.o sbexp.o sbfn.o sbfn1.o sbout.o \
sbraad.o sbrd1.o sbreed.o sbset.o sbtrp.o scale.o scraad.o scrd1.o \
screed.o sctrp.o sdeldp.o sdeps.o sdpsi.o seqchk.o shpcnt.o shpdim.o \
shppar.o shppr1.o sidtim.o sitcor.o sitpun.o sitred.o skycor.o skyin.o \
skypar.o skypun.o smlrot.o snitl.o snorm4.o snorml.o solctl.o solfrm.o \
solgr.o solprb.o soluse.o sotrp.o spotcd.o sptpun.o sptred.o staton.o \
stpctl.o subscp.o suicid.o syminv.o syminx.o tcheck.o terpf.o timdat.o \
timinc.o timrit.o timset.o tmend.o trgcnt.o trglic.o trgpic.o trplst.o \
trpnit.o uscosu.o ut1red.o ut2ut1.o uvectr.o vlrtrd.o wobble.o \
wobblf.o wobils.o wobred.o wupdat.o xa1ut1.o xsqrt.o xwobbl.o xxrd1.o \
xxrdbd.o ycoff.o ycoft.o yhermt.o ynutlb.o yplpcd.o yprtcd.o zenang.o \
zfill.o
pepsimp: $(SIMPOBJ) simpdum.ff
	cp simpdum.ff simpdum.f
	$(FORT) $(FFLAGS) -o pepsimp $(SIMPOBJ) simpdum.f
	rm simpdum.f

# include dependencies:
a1ut1.o: fcntrl.inc inodta.inc
a1ut1f.o: comdat.inc dtparm.inc empcnd.inc fcntrl.inc funcon.inc globdefs.inc leon.inc
acheck.o: fcntrl.inc
acmin.o: aprtbf.inc crdbuf.inc fcntrl.inc inodta.inc matzro.inc numnum.inc
acwrap.o: aprtbf.inc fcntrl.inc inodta.inc
adam.o: adams.inc coeff.inc globdefs.inc incon.inc inodta.inc nmstrt.inc output.inc stint.inc
addcln.o: fcntrl.inc filtim.inc inodta.inc
addinv.o: fcntrl.inc filtim.inc inodta.inc
adjast.o: adjstf.inc estvec.inc fcntrl.inc globdefs.inc rtsidesl.inc
adjbdy.o: adjstf.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc
adjfor.o: adjstf.inc globdefs.inc inodta.inc scoef4.inc
adjgrd.o: adjstf.inc globdefs.inc inodta.inc scoef4.inc
adjhar.o: adjstf.inc fcntrl.inc globdefs.inc inodta.inc
adjphs.o: adjstf.inc globdefs.inc inodta.inc phase.inc
adjrbs.o: adjstf.inc globdefs.inc inodta.inc rdbias.inc
adjsky.o: adjstf.inc globdefs.inc inodta.inc skystf.inc
adjspt.o: adjstf.inc funcon.inc globdefs.inc inodta.inc param.inc sptcrd.inc
adjust.o: adjstf.inc aprtbf.inc cureph.inc dtparm.inc empcnd.inc eqenox.inc ethhar.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc lcntrl.inc monhar.inc namtim.inc param.inc plnhar.inc psrstf.inc rtsidesl.inc scail.inc scoef4.inc stcord.inc
alaliz.o: fcntrl.inc obsdta.inc
albedo.o: globdefs.inc inodta.inc petuna.inc
analix.o: aprtbf.inc ciptrx.inc fcntrl.inc globdefs.inc inodta.inc mcnfrm.inc nrmmat.inc rtside.inc
analiz.o: aprtbf.inc ciptrx.inc fcntrl.inc inodta.inc nrmmat.inc obsdta.inc scail.inc
anaznm.o: inodta.inc
ancor.o: coord.inc fcntrl.inc globdefs.inc inodta.inc nutprc.inc param.inc sitcrd.inc
angaz.o: comdateq.inc coord.inc coordoeq.inc coordxeq.inc eqnphs.inc fcntrl.inc funcon.inc globdefs.inc jdfrequv.inc ltrapobs.inc number.inc nutprc.inc obscrd.inc param.inc sbcom.inc sitcrd.inc yvect.inc
angctl.o: comdateq.inc coord.inc coordoeq.inc coordxeq.inc fcntrl.inc globdefs.inc jdfrequv.inc kobequiv.inc nutprc.inc obscrd.inc param.inc sitcrd.inc spqind.inc yvect.inc
angdif.o: comdateq.inc coord.inc coordoeq.inc coordxeq.inc eqnphs.inc fcntrl.inc funcon.inc globdefs.inc jdfrequv.inc ltrapobs.inc number.inc nutprc.inc obscrd.inc param.inc sbcom.inc sitcrd.inc yvect.inc
angeq.o: comdateq.inc coord.inc coordoeq.inc coordxeq.inc eqnphs.inc fcntrl.inc funcon.inc globdefs.inc jdfrequv.inc ltrapobs.inc mnsprt.inc number.inc nutprc.inc obscrd.inc param.inc sbcom.inc yvect.inc
angmer.o: comdateq.inc coord.inc coordoeq.inc coordxeq.inc eqnphs.inc fcntrl.inc funcon.inc globdefs.inc jdfrequv.inc ltrapobs.inc number.inc nutprc.inc obscrd.inc param.inc sbcom.inc sitcrd.inc yvect.inc
angot.o: coord.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc loadid.inc mnsprt.inc namtim.inc number.inc nutprc.inc obscrd.inc param.inc redobs.inc sitcrd.inc yvect.inc
angout.o: coord.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc namtim.inc nutprc.inc obscrd.inc param.inc sitcrd.inc yvect.inc
angvp.o: comdateq.inc coord.inc coordxeq.inc globdefs.inc jdfrequv.inc number.inc obscrd.inc param.inc sbcom.inc sitcrd.inc
anset.o: anctrl.inc dtparm.inc eqenox.inc ethhar.inc fcntrl.inc globdefs.inc inodta.inc lcntrl.inc monhar.inc namtim.inc phase.inc plnhar.inc psrstf.inc rdbias.inc scoef4.inc skystf.inc sptcrd.inc stcord.inc
argmnt.o: funcon.inc
arrwrt.o: aprtbf.inc nrmmat.inc restor.inc
asbfnc.o: astroi.inc globdefs.inc intstf.inc petuna.inc sbstuf.inc sbthng.inc
asbset.o: astroi.inc funcon.inc globdefs.inc param.inc petuna.inc
astcrd.o: b2dtaint.inc b2ydta.inc inodta.inc prtpin.inc scdtaint.inc
asterp.o: fmmips.inc prtpin.inc
astrd1.o: b2dtaint.inc fcntrl.inc inodta.inc
atmden.o: globdefs.inc inodta.inc lothrf.inc param.inc petuna.inc sbthng.inc
b2rd1.o: b2dta.inc comdat.inc fcntrl.inc globdefs.inc inodta.inc
b2reed.o: b2dta.inc b2ydta.inc comdat.inc fcntrl.inc globdefs.inc inodta.inc sbdta.inc scdta.inc tapdta.inc tapdtp.inc
bcheck.o: fcntrl.inc
bdrd1.o: bdctrl.inc bddta.inc comdat.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc lcntrl.inc namtim.inc param.inc
bdreed.o: b2ydta.inc bdctrl.inc bddta.inc bdydta.inc comdat.inc comnut.inc fcntrl.inc globdefs.inc inodta.inc number.inc sbdta.inc scdta.inc tapdta.inc tapdte.inc tapdtm.inc tapdtp.inc trpcom.inc
bdrest.o: bdctrl.inc fcntrl.inc globdefs.inc inodta.inc rstart.inc
bdyprm.o: globdefs.inc
bdyred.o: empcnd.inc ethhar.inc france.inc globdefs.inc inodta.inc lcntrl.inc monhar.inc namtim.inc param.inc plndta.inc plnhar.inc psrstf.inc scoef4.inc smlbdy.inc
bihut1.o: argfun.inc fcntrl.inc inodta.inc
bihwob.o: fcntrl.inc funcon.inc inodta.inc nutprc.inc yvect.inc
blkdata1.o: funcon.inc globdefs.inc loadid.inc namtim.inc zeroes.inc
blkdata2.o: globdefs.inc incon.inc stint.inc
blkdata3.o: aacoff.inc orbstf.inc
blkdata4.o: prmnms.inc
blkdata7.o: comdateq.inc globdefs.inc namobs.inc
blkdata8.o: matzro.inc numnum.inc
blkdata9.o: rotcom.inc
bodast.o: bdctrl.inc empcnd.inc fmmips.inc funcon.inc globdefs.inc inodta.inc intstf.inc namtimq.inc param.inc petina.inc smlbdy.inc smlstf.inc
bodfn.o: aprtbf.inc bdctrl.inc bddtaint.inc bodstf.inc ellips.inc funcon.inc globdefs.inc intstf.inc orblun.inc output.inc param.inc prtcod.inc
bodgr.o: bodstf.inc globdefs.inc intstf.inc param.inc
bodhed.o: bdctrl.inc bddtaint.inc bodstf.inc fcntrl.inc globdefs.inc inodta.inc intstf.inc loadid.inc rstart.inc stint.inc zeroes.inc
bodntl.o: fcntrlx.inc globdefs.inc
bodout.o: adams.inc bdctrl.inc bddtaint.inc bdydta.inc bodstf.inc cassin.inc dumdum.inc empcnd.inc fcntrl.inc funcon.inc globdefs.inc incon.inc inodta.inc intstf.inc nmstrt.inc output.inc stint.inc
bodpun.o: fcntrl.inc globdefs.inc inodta.inc
bodred.o: fcntrlx.inc globdefs.inc inodta.inc maxkidat.inc param.inc
bodset.o: adams.inc bdctrl.inc bddtaint.inc bodstf.inc ellips.inc empcnd.inc funcon.inc globdefs.inc incon.inc inodta.inc intstf.inc metuna.inc namtim.inc nmstrt.inc param.inc petina.inc plndta.inc rstart.inc stint.inc
calcof.o: adams.inc coeff.inc globdefs.inc
change.o: ellips.inc funcon.inc
check.o: dtparm.inc empcnd.inc eqenox.inc ethhar.inc fcntrl.inc globdefs.inc inodta.inc lcntrl.inc maxmt0dt.inc monhar.inc namtim.inc obsdta.inc phase.inc plnhar.inc psrstf.inc rdbias.inc scoef4.inc skystf.inc sptcrd.inc stcord.inc
chncnc.o: funcon.inc
chncne.o: funcon.inc globdefs.inc param.inc
chncnt.o: fcntrl.inc globdefs.inc inodta.inc param.inc prtcod.inc
chncor.o: nutprc.inc
cmpar1.o: comdat.inc dltflg.inc dtparm.inc empcnd.inc fcntrl.inc globdefs.inc inodta.inc kobequiv.inc lcntrl.inc leon.inc loadid.inc ltrapx.inc mnsprt.inc mtrapx.inc namtim.inc number.inc obscrd.inc obsdta.inc obstap.inc param.inc sitcrd.inc skymap.inc spqind.inc stats.inc zeroes.inc
cmpar2.o: bdctrl.inc comdat.inc eqenox.inc eqnphs.inc fcntrl.inc globdefs.inc inodta.inc ltrapx.inc mnsprt.inc mtrapx.inc namtim.inc number.inc obscrd.inc obstap.inc phase.inc plndta.inc psrstf.inc rdbias.inc rotdta.inc sbdta.inc scdta.inc sitcrd.inc skymap.inc skystf.inc sptcrd.inc stats.inc stcord.inc tapdta.inc zeroes.inc
cmpar3.o: bdctrl.inc comdateq.inc coord.inc empcnd.inc eqnphs.inc fcntrl.inc globdefs.inc inodta.inc kobequiv.inc lcntrl.inc ltrapx.inc mnsprt.inc mtrapx.inc namtim.inc number.inc nutprc.inc obscrd.inc obstap.inc pemctl.inc plndta.inc rotdta.inc rtrdvl.inc sbdta.inc sitcrd.inc skymap.inc spqind.inc stats.inc tapdta.inc zeroes.inc
cmprst.o: comdat.inc fcntrl.inc globdefs.inc inodta.inc obsdta.inc zeroes.inc
cnstrn.o: fcntrl.inc
compar.o: bdctrl.inc comdat.inc crdbuf.inc fcntrl.inc globdefs.inc inodta.inc namtim.inc number.inc obsdta.inc obstap.inc plndta.inc rotdta.inc sbdta.inc scdta.inc sitcrd.inc stats.inc tapdta.inc zeroes.inc
comrit.o: comdat.inc eqnphs.inc fcntrl.inc globdefs.inc inodta.inc ltrapx.inc nutprc.inc obscrd.inc prpgat.inc sitcrd.inc zeroes.inc
comset.o: bdctrl.inc comdateq.inc empcnd.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc leon.inc ltrapx.inc namtim.inc number.inc obscrd.inc param.inc plndta.inc stats.inc timstf.inc
corchn.o: nutprc.inc
cpartl.o: comdateq.inc coord.inc eqnphs.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc kobequiv.inc ltrapobs.inc mnsprt.inc number.inc obscrd.inc param.inc partcm.inc pqind.inc radcrd.inc rotdta.inc rtrdvl.inc sbdta.inc sbdtavtp.inc scdta.inc scdtavtp.inc sitcrd.inc spqind.inc tabval.inc tapdte.inc tapdtm.inc tapdtp.inc trnocc.inc trpcom.inc
ctatf.o: comdateq.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc param.inc sitcrd.inc yvect.inc
ctreed.o: fcntrl.inc inodta.inc obscrd.inc tabval.inc
ctutf.o: dtparm.inc empcnd.inc fcntrl.inc funcon.inc globdefs.inc leon.inc
dcrfrm.o: fcntrl.inc globdefs.inc iptrst.inc mcnfrm.inc restor.inc rtsidesl.inc
decorr.o: bernum.inc fcntrl.inc globdefs.inc inodta.inc iptrst.inc ktrap.inc namobs.inc nrmmat.inc numnum.inc obsdta.inc restor.inc rtsidesl.inc scail.inc wrkcompr.inc
dguess.o: comdat.inc fcntrl.inc globdefs.inc inodta.inc obscrd.inc spqind.inc yvect.inc
difbdy.o: restor.inc
difeqn.o: anctrl.inc eqenox.inc eqenxm.inc globdefs.inc restor.inc wrkcomrs.inc
difhar.o: restor.inc
difnom.o: anctrl.inc dtparm.inc empcnd.inc ethhar.inc fcntrl.inc globdefs.inc inodta.inc lcntrl.inc mcnfrm.inc mhrfrm.inc monhar.inc namtim.inc param.inc plnhar.inc psrstf.inc psrstm.inc restor.inc scoef4.inc wrkcomrs.inc
difphs.o: anctrl.inc globdefs.inc phase.inc phasem.inc restor.inc wrkcomrs.inc
difpm4.o: restor.inc
difprm.o: restor.inc
difprt.o: fcntrl.inc inodta.inc
difrbs.o: anctrl.inc globdefs.inc rbiasm.inc rdbias.inc restor.inc wrkcomrs.inc
difsit.o: anctrl.inc globdefs.inc restor.inc stcord.inc stcrdm.inc wrkcomrs.inc
difspt.o: anctrl.inc globdefs.inc restor.inc sptcdm.inc sptcrd.inc wrkcomrs.inc
difstr.o: anctrl.inc globdefs.inc restor.inc skystf.inc skystm.inc wrkcomrs.inc
dish.o: dish1.inc
dishfc.o: dish1.inc
diurn.o: sitcrd.inc
diurnp.o: funcon.inc nutprc.inc
dlibra.o: funcon.inc
dltred.o: dltflg.inc fcntrl.inc inodta.inc obsdta.inc
dnpole.o: funcon.inc
dtpun.o: dtparm.inc inodta.inc
eatctl.o: coord.inc eqnphs.inc fcntrl.inc globdefs.inc number.inc param.inc prpgat.inc radcrd.inc
eatmmp.o: globdefs.inc obscrd.inc param.inc sitcrd.inc
eatzdl.o: fcntrl.inc globdefs.inc inodta.inc number.inc obscrd.inc param.inc sitcrd.inc
eatzmp.o: globdefs.inc param.inc sitcrd.inc
eclprc.o: aacoff.inc funcon.inc orbstf.inc precmn.inc
efgout.o: comdat.inc fcntrl.inc globdefs.inc inodta.inc nutprc.inc obscrd.inc param.inc sitcrd.inc yvect.inc
eioctl.o: coord.inc eqnphs.inc fcntrl.inc globdefs.inc inodta.inc number.inc obscrd.inc param.inc prpgat.inc radcrd.inc yvect.inc
elipt.o: ellips.inc
emipt.o: emcke.inc emmips.inc
emprd1.o: bddtaint.inc emmips.inc fcntrl.inc funcon.inc globdefs.inc intstf.inc metuna.inc morstf.inc namtim.inc param.inc plndta.inc tapdtplp.inc xprcom.inc yvectplp.inc yvectrd1.inc
emrd1.o: comdateq.inc empcnd.inc funcon.inc globdefs.inc namtim.inc param.inc pemctl.inc tapdta.inc trpcom.inc yvectrd1.inc
emreed.o: comdat.inc fcntrl.inc globdefs.inc inodta.inc number.inc plndta.inc tapdta.inc tapdte.inc
emtrp.o: comdateq.inc coord.inc fcntrl.inc globdefs.inc inodta.inc param.inc pqind.inc sitcrd.inc tapdta.inc tapdte.inc tapdtm.inc yvect.inc
eqnpun.o: eqenox.inc globdefs.inc inodta.inc
eqnred.o: eqenox.inc globdefs.inc inodta.inc
eqpcor.o: comdat.inc coord.inc eqnphs.inc fcntrl.inc funcon.inc globdefs.inc ltrapobs.inc number.inc obscrd.inc sitcrd.inc
errtot.o: comdat.inc fcntrl.inc globdefs.inc inodta.inc obscrd.inc obstuf.inc statsrad.inc timstf.inc
ertorb.o: cnthar.inc engstf.inc gleaks.inc globdefs.inc intstf.inc lothrf.inc output.inc petuna.inc sbstuf.inc sbthng.inc
ertpar.o: comdateq.inc coord.inc dtparm.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc leon.inc ltrapx.inc mtrapx.inc number.inc nutprc.inc obscrd.inc param.inc partcm.inc sitcrd.inc
ertrel.o: globdefs.inc intstf.inc sbstuf.inc sbthng.inc
ertsnt.o: bdctrl.inc bddta.inc bdydta.inc comdateq.inc globdefs.inc tabval.inc
eshape.o: comdateq.inc empcnd.inc eqnphs.inc funcon.inc globdefs.inc mnsprt.inc number.inc param.inc sitcrd.inc
etide.o: comdateq.inc coord.inc empcnd.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc obscrd.inc param.inc sitcrd.inc tidal.inc
eval.o: dumdum.inc globdefs.inc incon.inc inodta.inc output.inc stbtst.inc
evtrp.o: fcntrl.inc inodta.inc tabval.inc trpcom.inc
fcoher.o: filtim.inc inodta.inc
fercnt.o: coord.inc difnct.inc empcnd.inc fcntrl.inc globdefs.inc ltrapobs.inc number.inc obscrd.inc
ferctl.o: comdateq.inc coord.inc difnct.inc eqnphs.inc fcntrl.inc globdefs.inc kobequiv.inc ltrapobs.inc number.inc nutprc.inc obscrd.inc param.inc prpgat.inc scdta.inc sitcrd.inc yvect.inc
ferfrq.o: difnct.inc empcnd.inc fcntrl.inc globdefs.inc inodta.inc number.inc obscrd.inc sitcrd.inc
fermn.o: comdateq.inc coord.inc difnct.inc fcntrl.inc globdefs.inc mnsprt.inc namtim.inc number.inc nutprc.inc obscrd.inc param.inc radcrd.inc sbcom.inc sitcrd.inc yvect.inc
fermtr.o: comdat.inc difnct.inc fcntrl.inc globdefs.inc inodta.inc kobequiv.inc ltrapx.inc number.inc nutprc.inc obscrd.inc obstap.inc plndta.inc prpgat.inc redobs.inc sitcrd.inc spqind.inc statsrad.inc tidal.inc
fersb.o: comdateq.inc coord.inc difnct.inc fcntrl.inc globdefs.inc kobequiv.inc ltrapobs.inc mnsprt.inc number.inc nutprc.inc obscrd.inc param.inc radcrd.inc sbcom.inc scdta.inc sitcrd.inc yvect.inc
ferstr.o: comdateq.inc coord.inc difnct.inc fcntrl.inc globdefs.inc kobequiv.inc mnsprt.inc nutprc.inc obscrd.inc param.inc radcrd.inc sitcrd.inc yvect.inc
filctl.o: aprtbf.inc fcntrl.inc filptr.inc filtda.inc filtds.inc filtim.inc rtside.inc
filhed.o: filtds.inc inodta.inc rtside.inc
filstr.o: cureph.inc filtds.inc filtim.inc inodta.inc
filtep.o: fcntrl.inc filnit.inc filtds.inc filtim.inc globdefs.inc inodta.inc param.inc
filtin.o: filnit.inc filtda.inc filtds.inc filtim.inc inodta.inc
filtpn.o: filnit.inc filtim.inc inodta.inc
flatc.o: comdat.inc funcon.inc globdefs.inc shpcom.inc
flats.o: shpcom.inc
flured.o: fcntrl.inc globdefs.inc inodta.inc obscrd.inc param.inc sitcrd.inc
form.o: fcntrl.inc filptr.inc filstf.inc filtda.inc filtds.inc filtim.inc inodta.inc rtside.inc
formu.o: aprtbf.inc estvec.inc fcntrl.inc inodta.inc
frmbdy.o: restor.inc
frmeqn.o: anctrl.inc eqenox.inc eqenxm.inc globdefs.inc restor.inc wrkcomrs.inc
frmhar.o: restor.inc
frmhed.o: dtparm.inc empcnd.inc eqenxm.inc fcntrl.inc globdefs.inc inodta.inc lcntrl.inc mcnfrm.inc mhrfrm.inc namtim.inc phasem.inc psrstm.inc rbiasm.inc skystm.inc sptcdm.inc stcrdm.inc wrkcomrs.inc
frmmve.o: anctrl.inc dtparm.inc ethhar.inc globdefs.inc lcntrl.inc mcnfrm.inc mhrfrm.inc monhar.inc namtim.inc plnhar.inc psrstf.inc psrstm.inc restor.inc scoef4.inc zeroes.inc
frmphs.o: anctrl.inc globdefs.inc phase.inc phasem.inc restor.inc
frmprm.o: restor.inc
frmrbs.o: anctrl.inc globdefs.inc rbiasm.inc rdbias.inc restor.inc
frmsit.o: anctrl.inc globdefs.inc restor.inc stcord.inc stcrdm.inc
frmspt.o: anctrl.inc globdefs.inc restor.inc sptcdm.inc sptcrd.inc
frmstr.o: anctrl.inc globdefs.inc restor.inc skystf.inc skystm.inc
frstat.o: adjstf.inc fcntrl.inc globdefs.inc inodta.inc
funarg.o: argfun.inc
gausei.o: inodta.inc
getsol.o: fcntrl.inc filtda.inc filtim.inc inodta.inc rtside.inc
gpmpo.o: fcntrl.inc inodta.inc
grdpun.o: globdefs.inc inodta.inc
grdshp.o: comdat.inc globdefs.inc inodta.inc ltrapx.inc mtrapx.inc param.inc shpcom.inc shphar.inc
harshp.o: empcnd.inc funcon.inc globdefs.inc number.inc shpcom.inc shphar.inc
hdpni2.o: fcntrl.inc inodta.inc
hocred.o: comdateq.inc fcntrl.inc globdefs.inc hoccom.inc inodta.inc
horizn.o: fcntrl.inc globdefs.inc param.inc redobs.inc statsrad.inc
hpart2.o: globdefs.inc number.inc partcm.inc pemctl.inc tapdta.inc
hpartl.o: globdefs.inc maxkidat.inc partcm.inc pemctl.inc plndta.inc rotdta.inc sbdta.inc tapdta.inc
hrtgcn.o: cnthar.inc emmips.inc fcntrl.inc fmmips.inc globdefs.inc intstf.inc param.inc petuna.inc sbembr.inc sbrot.inc sbstuf.inc sbthng.inc trghar.inc yvectplp.inc
hrtrp.o: fcntrl.inc inodta.inc trpcom.inc
icnstr.o: filptr.inc filtda.inc
ierstide.o: funcon.inc globdefs.inc param.inc
imitl.o: emcke.inc emmips.inc
infred.o: globdefs.inc inodta.inc petuna.inc
inital.o: fcntrl.inc
initl.o: ellips.inc
input.o: crdbuf.inc fcntrl.inc inodta.inc nmlchrdt.inc
inwrit.o: france.inc globdefs.inc namtim.inc plndta.inc
jbdtst.o: b2dta.inc bdctrl.inc bddta.inc fcntrl.inc globdefs.inc inodta.inc namtim.inc plndta.inc
jnitl.o: funcon.inc
koutxpr.o: inodta.inc xprcom.inc
kp2ki.o: maxkidat.inc
legnds.o: globdefs.inc
lesin.o: funcon.inc inodta.inc
lesthp.o: funcon.inc globdefs.inc incon.inc intstf.inc param.inc petuna.inc sbstuf.inc sbthng.inc sscon.inc
lestin.o: thratt.inc
libchk.o: empcnd.inc globdefs.inc lcntrl.inc mtrapx.inc namtim.inc
linchk.o: adjstf.inc fcntrl.inc globdefs.inc inodta.inc
lmbcor.o: comdateq.inc coord.inc empcnd.inc funcon.inc globdefs.inc kobequiv.inc ltrapx.inc nutprc.inc obscrd.inc param.inc sitcrd.inc
lunorb.o: aacoff.inc funcon.inc orblun.inc orbstf.inc precmn.inc
lunset.o: aacoff.inc orbstf.inc precmn.inc
lvtbdy.o: comdat.inc globdefs.inc
lvthar.o: comdat.inc globdefs.inc
lvtprm.o: comdat.inc globdefs.inc
match.o: inodta.inc
media.o: comdateq.inc coord.inc funcon.inc globdefs.inc kobequiv.inc ltrapx.inc obscrd.inc param.inc sitcrd.inc yvect.inc
menmon.o: aacoff.inc fcntrl.inc funcon.inc
mensun.o: nutprc.inc
mnlib.o: comnut.inc
mnpcrd.o: aprtbf.inc bddtaint.inc bdydta.inc funcon.inc globdefs.inc inodta.inc param.inc prtcod.inc prtpin.inc scdtaint.inc
mnprd1.o: aprtbf.inc bddtaint.inc empcnd.inc fcntrl.inc globdefs.inc petina.inc prtpin.inc yvectrd1.inc
mnrd1.o: comdateq.inc comnut.inc globdefs.inc pemctl.inc tapdta.inc trpcom.inc yvectrd1.inc
mnreed.o: comdat.inc comnut.inc fcntrl.inc globdefs.inc inodta.inc number.inc plndta.inc tapdta.inc tapdtm.inc trpcom.inc
mnspt.o: comdateq.inc coord.inc empcnd.inc fcntrl.inc funcon.inc globdefs.inc hoccom.inc mnsprt.inc param.inc rotdtamr.inc
mntrp.o: comdateq.inc coord.inc empcnd.inc globdefs.inc param.inc pqind.inc tapdta.inc tapdtm.inc yvect.inc
monlib.o: angles.inc empcnd.inc funcon.inc globdefs.inc mnrtlb.inc precmn.inc prtcod.inc
monorb.o: empcnd.inc globdefs.inc sbthng.inc
monrel.o: bddtaint.inc globdefs.inc intstf.inc metuna.inc morcrd.inc morstf.inc param.inc prtcod.inc tapdtplp.inc xprcom.inc yvectplp.inc
monrot.o: cassin.inc empcnd.inc funcon.inc globdefs.inc inodta.inc mnrtlb.inc precmn.inc
montid.o: empcnd.inc globdefs.inc intstf.inc morcrd.inc morstf.inc nutces.inc
moon.o: aacoff.inc bdctrl.inc fcntrl.inc globdefs.inc inodta.inc metuna.inc plndta.inc
mored.o: empcnd.inc globdefs.inc harmor.inc intstf.inc mnrtlb.inc morcrd.inc morstf.inc param.inc
moredq.o: funcon.inc globdefs.inc harmor.inc morstf.inc param.inc
moredt.o: empcnd.inc globdefs.inc harmor.inc intstf.inc mnrtlb.inc morcrd.inc morstf.inc param.inc prtcod.inc
morffi.o: empcnd.inc globdefs.inc harmor.inc intstf.inc mnrtlb.inc morcrd.inc morstf.inc nutces.inc param.inc
morfn.o: ellips.inc emmips.inc empcnd.inc globdefs.inc harmor.inc inodta.inc intstf.inc metuna.inc mnrtlb.inc morcrd.inc morstf.inc orblun.inc output.inc param.inc precmn.inc prtcod.inc tapdtplp.inc xprcom.inc yvectplp.inc
morhar.o: globdefs.inc harmor.inc intstf.inc metuna.inc mnrtlb.inc morcrd.inc morstf.inc nutces.inc output.inc param.inc
morhed.o: empcnd.inc fcntrl.inc globdefs.inc inodta.inc intstf.inc loadid.inc metuna.inc morstf.inc namtim.inc param.inc plndta.inc rstart.inc stint.inc tapdtplp.inc zeroes.inc
morout.o: adams.inc bdydta.inc cassin.inc dumdum.inc ellips.inc empcnd.inc fcntrl.inc funcon.inc globdefs.inc incon.inc inodta.inc intstf.inc metuna.inc morstf.inc namtim.inc nmstrt.inc orblun.inc output.inc plndta.inc stint.inc zeroes.inc
morpar.o: empcnd.inc globdefs.inc harmor.inc intstf.inc morstf.inc param.inc
morset.o: adams.inc bddtaint.inc ellips.inc emmips.inc empcnd.inc ethhar.inc funcon.inc globdefs.inc harmor.inc incon.inc inodta.inc intstf.inc metuna.inc monhar.inc morstf.inc namtim.inc nmstrt.inc orblun.inc param.inc plndta.inc prtcod.inc rstart.inc stint.inc xprcom.inc
mrtpar.o: comdat.inc globdefs.inc ltrapx.inc mnsprt.inc mtrapx.inc number.inc partcm.inc pemctl.inc plndta.inc rotdtamr.inc tapdta.inc
mscfrc.o: mascn1.inc mascon.inc
mscset.o: funcon.inc globdefs.inc inodta.inc mascn1.inc mascon.inc param.inc
mtide.o: comdateq.inc coord.inc empcnd.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc mnsprt.inc obscrd.inc param.inc radcrd.inc tidal.inc
muler2.o: fcntrl.inc inodta.inc
mulerr.o: fcntrl.inc inodta.inc
mulitr.o: inodta.inc
multin.o: aprtbf.inc fcntrl.inc inodta.inc
mvmsrf.o: globdefs.inc petuna.inc
namprm.o: aprtbf.inc dtparm.inc empcnd.inc eqenox.inc ethhar.inc globdefs.inc inodta.inc lcntrl.inc monhar.inc namtim.inc param.inc plnhar.inc prmnms.inc psrstf.inc scoef4.inc skystf.inc stcord.inc
namred.o: globdefs.inc namtim.inc prmnms.inc
newhed.o: fcntrl.inc globdefs.inc inodta.inc mcnfrm.inc
newpg.o: fcntrl.inc inodta.inc
ninstp.o: adams.inc dumdum.inc globdefs.inc incon.inc inodta.inc nmstrt.inc output.inc stint.inc
nmlfix.o: nmlchrdt.inc
nomtrs.o: fcntrl.inc filtds.inc filtim.inc
norad.o: crdbuf.inc ebloc.inc funcon.inc globdefs.inc inodta.inc param.inc
nrmadd.o: fcntrl.inc inodta.inc
nrmcof.o: bernum.inc rtside.inc
nrmfrm.o: aprtbf.inc ciptrx.inc dltflg.inc dpsnec.inc fcntrl.inc globdefs.inc inodta.inc iptrst.inc mcnfrm.inc nrmwgt.inc obsdta.inc restor.inc rtside.inc
nrmict.o: aprtbf.inc bernum.inc empcnd.inc fcntrl.inc filmis.inc filtim.inc globdefs.inc inodta.inc ktrap.inc lcntrl.inc namobs.inc namtim.inc nrmmat.inc nrmwgt.inc obsdta.inc prdmov.inc rtside.inc wrkcompr.inc zeroes.inc
nrmrit.o: inodta.inc
nrmsav.o: aprtbf.inc cureph.inc fcntrl.inc filtim.inc globdefs.inc inodta.inc mcnfrm.inc nrmmat.inc restor.inc rtside.inc zeroes.inc
nrmscl.o: inodta.inc
nrmscq.o: inodta.inc
nrmset.o: fcntrl.inc nrmmat.inc numnum.inc rtside.inc
numint.o: adams.inc fcntrl.inc globdefs.inc incon.inc inodta.inc nmstrt.inc output.inc stint.inc
nuterp.o: prtpin.inc
nutrm.o: argfun.inc
obsbeg.o: comdat.inc eqnphs.inc fcntrl.inc globdefs.inc inodta.inc ltrapx.inc namtim.inc number.inc nutprc.inc obscrd.inc obstap.inc obstuf.inc redobs.inc sitcrd.inc skymap.inc statsrad.inc
obsrd1.o: dtparm.inc fcntrl.inc globdefs.inc inodta.inc mtrapx.inc obstap.inc
obsred.o: comdateq.inc empcnd.inc eqnphs.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc kobequiv.inc ltrapobs.inc mtrapx.inc namtim.inc number.inc obscrd.inc obsdta.inc obstap.inc param.inc prpgat.inc redobs.inc sitcrd.inc spqind.inc statsrad.inc trpcom.inc
occult.o: comdateq.inc coord.inc fcntrl.inc funcon.inc globdefs.inc jdfrequv.inc ltrapobs.inc mnsprt.inc namtim.inc number.inc nutprc.inc obscrd.inc param.inc sbcom.inc sitcrd.inc spqind.inc trnocc.inc watstf.inc yvect.inc
optic.o: comdateq.inc coord.inc empcnd.inc eqnphs.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc ltrapobs.inc number.inc obscrd.inc obstap.inc param.inc redobs.inc sitcrd.inc statsrad.inc
partl.o: bddta.inc cmcke.inc comdateq.inc coord.inc emmips.inc eqnphs.inc fcntrl.inc funcon.inc globdefs.inc kobequiv.inc lfix.inc ltrap.inc mnsprt.inc mtrap.inc namtim.inc number.inc obscrd.inc param.inc partcm.inc pemctl.inc radcrd.inc rotdtamr.inc rtrdvl.inc sbcom.inc sbdta.inc scdta.inc sitcrd.inc tapdta.inc tidal.inc trnocc.inc yvect.inc
partl1.o: atmos.inc comdateq.inc coord.inc difnct.inc emmips.inc empcnd.inc eqnphs.inc fcntrl.inc funcon.inc globdefs.inc kobequiv.inc ltrap.inc mnsprt.inc mtrap.inc namtim.inc number.inc obscrd.inc param.inc partcm.inc pemctl.inc radcrd.inc sbdta.inc scdta.inc tapdta.inc trnocc.inc yvect.inc
partvl.o: fcntrl.inc partcm.inc rtrdvl.inc
patctl.o: atmos.inc comdat.inc coord.inc empcnd.inc fcntrl.inc globdefs.inc inodta.inc mnsprt.inc number.inc obscrd.inc param.inc radcrd.inc scdta.inc
pbdic.o: globdefs.inc inodta.inc namtim.inc
pcheck.o: bdctrl.inc fcntrl.inc france.inc globdefs.inc inodta.inc lcntrl.inc namtim.inc
pcops.o: globdefs.inc lfix.inc ltrapobs.inc mtrap.inc partcm.inc
pepmain.o: bdctrl.inc crdbuf.inc fcntrl.inc globdefs.inc inodta.inc namtim.inc timstf.inc
pepnut.o: funcon.inc
peptic.o: crdbuf.inc nmlchrdt.inc
peptin.o: crdbuf.inc
perchk.o: aprtbf.inc bdctrl.inc fcntrl.inc filtda.inc filtds.inc france.inc globdefs.inc inodta.inc obsdta.inc plndta.inc
perdta.o: aprtbf.inc globdefs.inc inodta.inc maxmt0dt.inc obsdta.inc plndta.inc
pexlic.o: comdat.inc empcnd.inc ennips.inc funcon.inc globdefs.inc lcntrl.inc ltrapx.inc maxpxdat.inc mtrapx.inc namtim.inc number.inc obscrd.inc zeroes.inc
phadop.o: comdat.inc coord.inc difnct.inc globdefs.inc kobequiv.inc ltrapobs.inc number.inc obscrd.inc param.inc sitcrd.inc spqind.inc yvect.inc
phscor.o: comdateq.inc coord.inc globdefs.inc partcm.inc trnocc.inc
phsprm.o: globdefs.inc phase.inc
phspun.o: globdefs.inc inodta.inc phase.inc
phsred.o: globdefs.inc inodta.inc namtim.inc phase.inc
planet.o: aacoff.inc aprtbf.inc bdctrl.inc empcnd.inc fcntrl.inc globdefs.inc inodta.inc namtimq.inc param.inc petina.inc petuna.inc plndta.inc stbtst.inc stint.inc timstf.inc
plcmc.o: b2dta.inc b2ydta.inc bddta.inc bdydta.inc cmcke.inc comdat.inc fcntrl.inc globdefs.inc inodta.inc namtim.inc number.inc param.inc tabval.inc
plcms.o: b2dta.inc bddta.inc cmcke.inc comdat.inc empcnd.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc maxpcmdt.inc namtim.inc number.inc param.inc plndta.inc stats.inc
plink.o: inodta.inc loadid.inc
plnhed.o: fcntrl.inc globdefs.inc inodta.inc intstf.inc loadid.inc param.inc petuna.inc rstart.inc sbstuf.inc stint.inc tapdtplp.inc zeroes.inc
plnorb.o: astroi.inc cnthar.inc drgprt.inc engstf.inc gleaks.inc globdefs.inc intstf.inc lothrf.inc param.inc petuna.inc sbstuf.inc sbthng.inc
plnsat.o: fcntrl.inc globdefs.inc intstf.inc param.inc petuna.inc sbstuf.inc sbthng.inc
plpcrd.o: emcke.inc emmips.inc funcon.inc globdefs.inc inodta.inc orblun.inc prtcod.inc prtpin.inc tapdtplp.inc yvectplp.inc
plpndx.o: globdefs.inc inodta.inc
plprd1.o: bddtaint.inc fcntrl.inc globdefs.inc intstf.inc namtimq.inc petuna.inc plndta.inc prtpin.inc sbstuf.inc tapdtplp.inc xprcom.inc yvectplp.inc yvectrd1.inc
plptrp.o: globdefs.inc prtpin.inc tapdtplp.inc yvectplp.inc
plrd1.o: comdateq.inc empcnd.inc fcntrl.inc funcon.inc globdefs.inc lcntrl.inc namtim.inc number.inc param.inc pemctl.inc plnhar.inc rotcom.inc scoef4.inc shpcom.inc shphar.inc tapdta.inc trpcom.inc yvectrd1.inc
plreed.o: comdat.inc fcntrl.inc globdefs.inc inodta.inc namtim.inc number.inc tapdta.inc tapdtp.inc trpcom.inc
pltrp.o: comdateq.inc coord.inc emcke.inc emmips.inc fcntrl.inc globdefs.inc pemctl.inc pqind.inc tapdta.inc tapdtp.inc yvect.inc
pnitl.o: funcon.inc
pprctl.o: aprtbf.inc ciptrx.inc fcntrl.inc globdefs.inc inodta.inc iptrst.inc mcnfrm.inc nrmmat.inc numnum.inc restor.inc rtside.inc scail.inc zeroes.inc
pprfrm.o: restor.inc
prces.o: fcntrl.inc precmn.inc
prces9.o: funcon.inc
prcnt9.o: funcon.inc
prcnut.o: argfun.inc comdateq.inc comnut.inc empcnd.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc leon.inc nutprc.inc
prdbdy.o: wrkcompm.inc
prdeqs.o: anctrl.inc bernum.inc dtparm.inc eqenox.inc ethhar.inc fcntrl.inc globdefs.inc inodta.inc ktrap.inc lcntrl.inc monhar.inc phase.inc plnhar.inc psrstf.inc rdbias.inc scoef4.inc skystf.inc sptcrd.inc stcord.inc wrkcompm.inc
prdgrd.o: bernum.inc globdefs.inc ktrap.inc wrkcompm.inc
prdhar.o: wrkcompm.inc
prdict.o: bernum.inc fcntrl.inc globdefs.inc inodta.inc ktrap.inc obsdta.inc prdmov.inc rtsidesl.inc wrkcompr.inc
prdobs.o: bernum.inc dltflg.inc fcntrl.inc globdefs.inc inodta.inc ktrap.inc loadid.inc namobs.inc nrmwgt.inc obsdta.inc prdmov.inc wrkcompr.inc zeroes.inc
prdprm.o: wrkcompm.inc
prdset.o: bernum.inc eqenox.inc globdefs.inc inodta.inc ktrap.inc namtim.inc phase.inc psrstf.inc rdbias.inc skystf.inc sptcrd.inc stcord.inc wrkcompr.inc
prdsid.o: bernum.inc fcntrl.inc globdefs.inc ktrap.inc nrmmat.inc prdmov.inc rtsidesl.inc
prdvct.o: wrkcompm.inc
preces.o: empcnd.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc leon.inc mnsprt.inc number.inc nutprc.inc
prenut.o: nutces.inc precmn.inc prtcod.inc
prepda.o: fcntrl.inc filtda.inc globdefs.inc
prmred.o: aprtbf.inc bdctrl.inc dtparm.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc lcntrl.inc loadid.inc mascon.inc maxmt0dt.inc nrmwgt.inc obsdta.inc param.inc plndta.inc prpgat.inc timstf.inc
prmtrn.o: aprtbf.inc correl.inc fcntrl.inc inodta.inc maxcrrdt.inc
prnacm.o: aprtbf.inc fcntrl.inc inodta.inc
prncrd.o: fcntrl.inc globdefs.inc inodta.inc namtim.inc sptcrd.inc stcord.inc
prnfil.o: fcntrl.inc filnit.inc filtda.inc filtds.inc filtim.inc inodta.inc
prnobs.o: dltflg.inc eqenox.inc fcntrl.inc globdefs.inc inodta.inc namtim.inc obsdta.inc phase.inc rdbias.inc skystf.inc sptcrd.inc stcord.inc
prnpsr.o: fcntrl.inc globdefs.inc inodta.inc psrstf.inc sptcrd.inc
prnstr.o: fcntrl.inc globdefs.inc inodta.inc skystf.inc
prnsum.o: aprtbf.inc fcntrl.inc filstf.inc filtds.inc filtim.inc inodta.inc revptr.inc
prntot.o: aprtbf.inc bdctrl.inc dtparm.inc empcnd.inc ethhar.inc fcntrl.inc france.inc funcon.inc globdefs.inc inodta.inc lcntrl.inc monhar.inc namtim.inc nrmwgt.inc obsdta.inc param.inc plndta.inc plnhar.inc scoef4.inc smlbdy.inc timstf.inc
propco.o: fcntrl.inc inodta.inc kobequiv.inc obscrd.inc prpgat.inc
prpfrq.o: obscrd.inc sitcrd.inc spqind.inc
prplog.o: fcntrl.inc inodta.inc obscrd.inc prpgat.inc
prpstr.o: difnct.inc fcntrl.inc inodta.inc number.inc obscrd.inc prpgat.inc
prpzen.o: coord.inc fcntrl.inc funcon.inc inodta.inc number.inc obscrd.inc prpgat.inc radcrd.inc sitcrd.inc
prpzro.o: prpgat.inc
prtcrd.o: bddtaint.inc bdydta.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc orblun.inc param.inc prtcod.inc prtpin.inc scdtaint.inc
prterp.o: prtcod.inc prtpin.inc
prtpr1.o: ltrapobs.inc mnsprt.inc partcm.inc prtpr9.inc
prtrd1.o: bdctrl.inc bddtaint.inc fcntrl.inc globdefs.inc inodta.inc param.inc prtpin.inc
prtrel.o: comdateq.inc globdefs.inc ltrapobs.inc obscrd.inc partcm.inc
psrpar.o: comdateq.inc empcnd.inc ennips.inc eqnphs.inc funcon.inc globdefs.inc ltrapx.inc mtrapx.inc number.inc obscrd.inc param.inc prtpr9.inc psrcom.inc rtrdvl.inc yvect.inc
psrphs.o: comdateq.inc eqnphs.inc globdefs.inc obscrd.inc psrcom.inc
psrpun.o: globdefs.inc inodta.inc psrstf.inc
punch.o: empcnd.inc ethhar.inc fcntrl.inc globdefs.inc inodta.inc lcntrl.inc monhar.inc namtim.inc param.inc plndta.inc plnhar.inc scoef4.inc zeroes.inc
pvbus.o: comdateq.inc coord.inc empcnd.inc globdefs.inc obscrd.inc
pvcrd.o: comdateq.inc coord.inc empcnd.inc funcon.inc globdefs.inc mnsprt.inc number.inc param.inc sbcom.inc scdta.inc
qset.o: rotcom.inc
radar.o: comdateq.inc difnct.inc empcnd.inc eqnphs.inc fcntrl.inc globdefs.inc inodta.inc kobequiv.inc lfix.inc ltrapx.inc mtrapx.inc number.inc nutprc.inc obscrd.inc obstap.inc obstuf.inc param.inc plndta.inc prpgat.inc redobs.inc sitcrd.inc spqind.inc statsrad.inc tidal.inc
radctl.o: comdateq.inc coord.inc difnct.inc eqnphs.inc fcntrl.inc funcon.inc globdefs.inc kobequiv.inc ltrapobs.inc mnsprt.inc number.inc nutprc.inc obscrd.inc param.inc plnhar.inc prpgat.inc sitcrd.inc spqind.inc yvect.inc
radext.o: comdateq.inc coord.inc eqnphs.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc kobequiv.inc ltrapobs.inc mnsprt.inc number.inc nutprc.inc obscrd.inc param.inc plnhar.inc sitcrd.inc spqind.inc yvect.inc
radgps.o: comdateq.inc coord.inc eqnphs.inc fcntrl.inc globdefs.inc kobequiv.inc ltrapobs.inc mnsprt.inc number.inc nutprc.inc obscrd.inc param.inc plnhar.inc prpgat.inc sbcom.inc sitcrd.inc spqind.inc yvect.inc
radlnd.o: comdateq.inc coord.inc eqnphs.inc globdefs.inc kobequiv.inc ltrapobs.inc mnsprt.inc number.inc nutprc.inc obscrd.inc param.inc prpgat.inc sitcrd.inc spqind.inc yvect.inc
radmn.o: comdateq.inc coord.inc eqnphs.inc fcntrl.inc globdefs.inc kobequiv.inc ltrapobs.inc mnsprt.inc number.inc nutprc.inc obscrd.inc param.inc plnhar.inc prpgat.inc sbcom.inc sitcrd.inc spqind.inc yvect.inc
radnrm.o: comdateq.inc coord.inc eqnphs.inc globdefs.inc kobequiv.inc ltrapobs.inc number.inc obscrd.inc param.inc pqind.inc prpgat.inc sbcom.inc sitcrd.inc yvect.inc
radns.o: comdateq.inc coord.inc eqnphs.inc globdefs.inc kobequiv.inc ltrapobs.inc number.inc nutprc.inc obscrd.inc param.inc prpgat.inc sbcom.inc sitcrd.inc spqind.inc yvect.inc
radpl.o: comdateq.inc coord.inc eqnphs.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc kobequiv.inc ltrapobs.inc mnsprt.inc number.inc nutprc.inc obscrd.inc param.inc plnhar.inc prpgat.inc sitcrd.inc spqind.inc yvect.inc
radpsr.o: comdateq.inc coord.inc difnct.inc ennips.inc eqnphs.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc kobequiv.inc ltrapx.inc mnsprt.inc number.inc obscrd.inc param.inc prpgat.inc rtrdvl.inc sitcrd.inc yvect.inc
radrel.o: comdateq.inc coord.inc difnct.inc eqnphs.inc fcntrl.inc globdefs.inc inodta.inc kobequiv.inc ltrapobs.inc obscrd.inc param.inc sitcrd.inc yvect.inc
radsb.o: comdateq.inc coord.inc eqnphs.inc globdefs.inc inodta.inc kobequiv.inc ltrapobs.inc number.inc nutprc.inc obscrd.inc param.inc prpgat.inc sbcom.inc sitcrd.inc spqind.inc yvect.inc
radvaf.o: funcon.inc
rbsprm.o: globdefs.inc rdbias.inc
rbspun.o: globdefs.inc inodta.inc rdbias.inc
rbsred.o: globdefs.inc inodta.inc namtim.inc rdbias.inc
reflcc.o: funcon.inc inodta.inc
refrct.o: fcntrl.inc inodta.inc obscrd.inc sitcrd.inc
restrt.o: fcntrl.inc globdefs.inc inodta.inc rstart.inc
rfilda.o: filtda.inc
rfiltr.o: filtds.inc
rinsne.o: filtds.inc
rngmod.o: fcntrl.inc globdefs.inc inodta.inc ltrapx.inc obscrd.inc
rochng.o: funcon.inc
rotab.o: rotcom.inc
rotat.o: inodta.inc
rotcrd.o: globdefs.inc pqind.inc rotdta.inc
rotips.o: rotcom.inc
rotmat.o: cnthar.inc funcon.inc globdefs.inc intstf.inc mnrtlb.inc nutces.inc petuna.inc rtpars.inc sbstuf.inc trghar.inc
rotset.o: funcon.inc
rroad.o: adams.inc coeff.inc globdefs.inc incon.inc inodta.inc nmstrt.inc output.inc stint.inc
rsmat.o: filtda.inc filtds.inc filtim.inc inodta.inc
rtrd1.o: comdateq.inc globdefs.inc namtim.inc number.inc rotdta.inc trpcom.inc yvectrd1.inc
rtreed.o: globdefs.inc inodta.inc rotdta.inc trpcom.inc
savhed.o: dtparm.inc empcnd.inc eqenox.inc ethhar.inc fcntrl.inc globdefs.inc lcntrl.inc loadid.inc monhar.inc namtim.inc param.inc phase.inc plnhar.inc psrstf.inc rdbias.inc scoef4.inc skystf.inc sptcrd.inc stcord.inc zeroes.inc
sbeng.o: engstf.inc inodta.inc
sbexp.o: empcnd.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc intstf.inc namtimq.inc param.inc petuna.inc sbembr.inc sbrot.inc sbstuf.inc sbthng.inc
sbfn.o: cnthar.inc ellips.inc fcntrl.inc fmmips.inc globdefs.inc inodta.inc intstf.inc orblun.inc output.inc param.inc petuna.inc prtcod.inc sbembr.inc sbrot.inc sbstuf.inc sbthng.inc stint.inc tapdtplp.inc trghar.inc xprcom.inc yvectplp.inc
sbfn1.o: cnthar.inc drgprt.inc ellips.inc fcntrl.inc fmmips.inc globdefs.inc incon.inc inodta.inc intstf.inc mascn1.inc namtimq.inc param.inc petuna.inc prtcod.inc rtpars.inc sbembr.inc sbrot.inc sbstuf.inc sbthng.inc tapdtplp.inc trghar.inc xprcom.inc yvectplp.inc
sbfout.o: adams.inc dumdum.inc ellips.inc fcntrl.inc filtds.inc filtim.inc funcon.inc globdefs.inc incon.inc inodta.inc intstf.inc lothrf.inc nmstrt.inc output.inc param.inc petuna.inc rstart.inc sbstuf.inc sscon.inc stbtst.inc stint.inc zeroes.inc
sbout.o: adams.inc dumdum.inc ellips.inc fcntrl.inc globdefs.inc incon.inc inodta.inc intstf.inc nmstrt.inc output.inc param.inc petuna.inc rstart.inc sbstuf.inc sscon.inc stint.inc zeroes.inc
sbraad.o: fcntrl.inc flcomp.inc globdefs.inc inodta.inc number.inc sbdta.inc sbdtavtp.inc spqind.inc trpcom.inc
sbrd1.o: comdateq.inc emcke.inc empcnd.inc flcomp.inc funcon.inc globdefs.inc namtim.inc number.inc param.inc sbdta.inc trpcom.inc yvectrd1.inc
sbreed.o: fcntrl.inc globdefs.inc inodta.inc number.inc sbdta.inc trpcom.inc
sbrset.o: lothrf.inc
sbset.o: adams.inc astroi.inc b2dtaint.inc bddtaint.inc cnthar.inc drgprt.inc ellips.inc emmips.inc empcnd.inc ethhar.inc fcntrl.inc fmmips.inc funcon.inc gleaks.inc globdefs.inc incon.inc inodta.inc intstf.inc lothrf.inc monhar.inc namtimq.inc nmstrt.inc param.inc petina.inc petuna.inc plnhar.inc prtpin.inc rstart.inc rtpars.inc sbstuf.inc scoef4.inc smlbdy.inc smlstf.inc sscon.inc stint.inc thratt.inc trghar.inc xprcom.inc
sbtrp.o: comdateq.inc emcke.inc emmips.inc globdefs.inc sbcom.inc sbdta.inc sbdtavtp.inc spqind.inc
scale.o: funcon.inc
scraad.o: globdefs.inc inodta.inc scdta.inc scdtavtp.inc trpcom.inc
scrd1.o: globdefs.inc namtim.inc number.inc scdta.inc trpcom.inc yvectrd1.inc
screed.o: globdefs.inc inodta.inc namtim.inc number.inc scdta.inc trpcom.inc
sctrp.o: globdefs.inc sbcom.inc scdta.inc scdtavtp.inc spqind.inc trpcom.inc
seqchk.o: comdat.inc fcntrl.inc globdefs.inc inodta.inc namtim.inc obstap.inc stats.inc
sgp.o: ebloc.inc funcon.inc
shadow.o: funcon.inc globdefs.inc lothrf.inc param.inc prtcod.inc sbthng.inc
shpcnt.o: fcntrl.inc globdefs.inc ltrapx.inc monhar.inc mtrapx.inc number.inc shphar.inc
shppar.o: comdat.inc globdefs.inc ltrapx.inc mtrapx.inc shpcom.inc shphar.inc
shppr1.o: comdateq.inc globdefs.inc ltrapx.inc param.inc shpcom.inc shphar.inc
shpprm.o: globdefs.inc
sidtim.o: empcnd.inc fcntrl.inc funcon.inc globdefs.inc
sitcor.o: comdateq.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc nutprc.inc param.inc sitcrd.inc
sitpun.o: globdefs.inc inodta.inc stcord.inc
sitred.o: empcnd.inc funcon.inc globdefs.inc inodta.inc stcord.inc
sj2set.o: funcon.inc globdefs.inc intstf.inc param.inc
skale.o: funcon.inc globdefs.inc param.inc
skycor.o: funcon.inc globdefs.inc ltrapx.inc shpcom.inc skymap.inc
skyin.o: globdefs.inc inodta.inc skystf.inc
skypar.o: comdat.inc funcon.inc globdefs.inc ltrapx.inc mtrapx.inc shpcom.inc
skypun.o: globdefs.inc inodta.inc skystf.inc
snorm2.o: inodta.inc
snorm4.o: inodta.inc
snorm5.o: fcntrl.inc filtim.inc inodta.inc nrmmat.inc rtside.inc
snorm6.o: fcntrl.inc
snorml.o: aprtbf.inc correl.inc fcntrl.inc filtim.inc inodta.inc nrmmat.inc rtside.inc
solctl.o: aprtbf.inc ciptrx.inc fcntrl.inc inodta.inc obsdta.inc rtsidesl.inc scail.inc
solfrm.o: fcntrl.inc globdefs.inc inodta.inc mcnfrm.inc nrmmat.inc numnum.inc restor.inc rtside.inc
solgr.o: bddtaint.inc globdefs.inc intstf.inc param.inc petuna.inc prtcod.inc sbembr.inc sbstuf.inc sbthng.inc xprcom.inc
solppr.o: fcntrl.inc globdefs.inc inodta.inc iptrst.inc mcnfrm.inc nrmmat.inc restor.inc rtsidesl.inc
solprb.o: astroi.inc globdefs.inc intstf.inc param.inc petuna.inc sbrot.inc sbstuf.inc sbthng.inc smlbdy.inc smlstf.inc trghar.inc xprcom.inc
soluse.o: fcntrl.inc inodta.inc rtsidesl.inc scail.inc
sotrp.o: bdctrl.inc bddta.inc bdydta.inc cmcke.inc globdefs.inc number.inc param.inc trpcom.inc
spotcd.o: comdateq.inc eqnphs.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc mnsprt.inc number.inc nutprc.inc param.inc prtpr9.inc
spotpl.o: comdateq.inc empcnd.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc leon.inc mnsprt.inc number.inc obscrd.inc param.inc partcm.inc rotcom.inc
sptprm.o: globdefs.inc sptcrd.inc
sptpun.o: globdefs.inc inodta.inc sptcrd.inc
sptred.o: empcnd.inc funcon.inc globdefs.inc inodta.inc namtim.inc sptcrd.inc
stangl.o: comdateq.inc coord.inc eqnphs.inc fcntrl.inc funcon.inc globdefs.inc kobequiv.inc ltrapobs.inc mnsprt.inc number.inc obscrd.inc param.inc sbcom.inc sitcrd.inc spqind.inc yvect.inc
staton.o: fcntrl.inc funcon.inc sscon.inc
stdldp.o: comdateq.inc coord.inc eqnphs.inc globdefs.inc kobequiv.inc ltrapobs.inc mnsprt.inc number.inc obscrd.inc param.inc sbcom.inc sitcrd.inc spqind.inc yvect.inc
stmdld.o: comdateq.inc coord.inc eqnphs.inc globdefs.inc kobequiv.inc ltrapobs.inc namtim.inc number.inc nutprc.inc obscrd.inc param.inc sbcom.inc sitcrd.inc spqind.inc yvect.inc
stmrdr.o: comdat.inc fcntrl.inc globdefs.inc inodta.inc kobequiv.inc ltrapx.inc namtim.inc number.inc nutprc.inc obscrd.inc obstap.inc param.inc redobs.inc sitcrd.inc spqind.inc statsrad.inc
stoptc.o: comdat.inc coord.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc ltrapx.inc namtim.inc number.inc obscrd.inc obstap.inc param.inc redobs.inc sitcrd.inc statsrad.inc
stpctl.o: globdefs.inc stint.inc
stradr.o: fcntrl.inc globdefs.inc inodta.inc kobequiv.inc ltrapx.inc namtim.inc number.inc obscrd.inc obstap.inc param.inc redobs.inc spqind.inc statsrad.inc
subscp.o: comdat.inc coord.inc fcntrl.inc funcon.inc globdefs.inc mnsprt.inc number.inc obscrd.inc sbcom.inc yvect.inc
suicid.o: fcntrl.inc inodta.inc
sxy06.o: funcon.inc
tcheck.o: inodta.inc
timdat.o: timstf.inc
timrit.o: inodta.inc timstf.inc
timset.o: timstf.inc
tmend.o: fcntrl.inc inodta.inc timstf.inc
trgcnt.o: ethhar.inc fcntrl.inc globdefs.inc ltrapx.inc monhar.inc mtrapx.inc number.inc pemctl.inc plnhar.inc sbdta.inc scoef4.inc
trglic.o: comdat.inc globdefs.inc inodta.inc lcntrl.inc ltrapx.inc mtrapx.inc namtim.inc number.inc plnhar.inc scoef4.inc
trgpic.o: comdat.inc globdefs.inc ltrapx.inc mtrapx.inc number.inc partcm.inc pemctl.inc tapdta.inc
trnsit.o: comdateq.inc empcnd.inc fcntrl.inc funcon.inc globdefs.inc inodta.inc ltrapobs.inc namtim.inc number.inc obscrd.inc obstap.inc param.inc redobs.inc sitcrd.inc statsrad.inc stcord.inc trnocc.inc watstf.inc
trplst.o: fcntrl.inc inodta.inc
trpnit.o: bddta.inc trpcom.inc
uscosu.o: fcntrl.inc inodta.inc prpgat.inc
use.o: fcntrl.inc filstf.inc filtda.inc filtds.inc filtim.inc inodta.inc revptr.inc rtside.inc scail.inc
ut1red.o: argfun.inc empcnd.inc fcntrl.inc globdefs.inc inodta.inc leon.inc obscrd.inc
ut2ut1.o: funcon.inc
uvectr.o: fcntrl.inc inodta.inc
vlrtrd.o: coord.inc fcntrl.inc rtrdvl.inc
wfilda.o: filtda.inc
wfilde.o: filtda.inc
wfiltr.o: filtds.inc
winsne.o: aprtbf.inc fcntrl.inc filmis.inc filtds.inc filtim.inc inodta.inc
wobble.o: fcntrl.inc
wobblf.o: dtparm.inc fcntrl.inc leon.inc
wobred.o: fcntrl.inc inodta.inc obscrd.inc
wotsne.o: filtda.inc filtds.inc
xxrd1.o: empcnd.inc fcntrl.inc globdefs.inc inodta.inc lcntrl.inc namtim.inc yvectrd1.inc
xxrdbd.o: b2dta.inc bddta.inc comdat.inc empcnd.inc globdefs.inc namtim.inc plndta.inc yvectrd1.inc
ycoff.o: tabval.inc
ycoft.o: tabval.inc trpcom.inc
yprtcd.o: prtpin.inc
