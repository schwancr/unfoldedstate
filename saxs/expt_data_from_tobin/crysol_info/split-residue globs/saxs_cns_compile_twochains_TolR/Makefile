########################################################
#       Windows-NT   Tue Aug 31 11:59:33 EDT 1999      #
########################################################
    
OBJECTS = angledb.obj \
aria.obj \
ariass.obj \
arical.obj \
ariman.obj \
arivio.obj \
array.obj \
carb.obj \
cluster.obj \
cns.obj \
connect.obj \
coorio.obj \
corman.obj \
coup.obj \
cstack.obj \
cstran.obj \
dcartesian.obj \
define.obj \
diff_anis.obj \
distan.obj \
dsearch.obj \
dtorsion_md.obj \
dtorsion_top.obj \
dynio.obj \
eangle.obj \
ebond.obj \
enbond.obj \
energy.obj \
etor.obj \
expression.obj \
fastnb.obj \
fft.obj \
fftvfy.obj \
flagmap.obj \
genic.obj \
heap.obj \
heapvfy.obj \
hbdb.obj \
saxs.obj \
dcsa.obj \
collapse.obj \
imf.obj \
initia.obj \
lbfgs.obj \
lidens.obj \
machine_f.obj \
machvar.obj \
mapyard.obj \
matrix.obj \
mindist.obj \
mmdisg.obj \
modes.obj \
mrigid.obj \
mtfio.obj \
nbonds.obj \
ncs.obj \
noe.obj \
onebond.obj \
parmio.obj \
parser.obj \
pick.obj \
pig.obj \
plane.obj \
powell.obj \
print.obj \
prot.obj \
psearch.obj \
psfio.obj \
rama_gaussians.obj \
read_write.obj \
rotate.obj \
rotlsq.obj \
rotman.obj \
rsearch.obj \
rtfio.obj \
segmnt.obj \
selrpn.obj \
smf.obj \
sort.obj \
string.obj \
surfac.obj \
susc_anis.obj \
test.obj \
tsmap.obj \
tsmapd.obj \
tsmapdlm.obj \
tsmaplm.obj \
tsmaps.obj \
update.obj \
util.obj \
vector.obj \
vector_func.obj \
xasymm.obj \
xdeclare.obj \
xdo.obj \
xdofft.obj \
xdofried.obj \
xfalgebra.obj \
xfft.obj \
xfunct.obj \
xmalgebra.obj \
xmask.obj \
xmaskio.obj \
xmaverage.obj \
xmaxl.obj \
xmphisto.obj \
xmpsele.obj \
xmread.obj \
xmulti.obj \
xmwrite.obj \
xparse.obj \
xpeakpik.obj \
xpred.obj \
xprobab.obj \
xprox.obj \
xray.obj \
xreduce.obj \
xreflstar.obj \
xremap.obj \
xrmani.obj \
xropti.obj \
xscale.obj \
xsdev.obj \
xsfsele.obj \
xshow.obj \
xsigmaa.obj \
xskel.obj \
xspecial.obj \
xstats.obj \
xsymm.obj \
xtarget.obj \
xutil.obj \
xyzparse.obj \
machine_c.obj \
orientations.obj \
dmemory.obj
 
###############################################################################
#                                                                             #
#       Machine dependent Makefile header for CNSsolve                        #
#       Paul Adams  20-10-98                                                  #
#       copyright Yale University                                             #
#                                                                             #
###############################################################################

TMP=.

cns_solve.exe: $(OBJECTS)
	df /nologo /fpscomp /exe:cns_solve.exe $(OBJECTS)

.F.obj:
	df  -optimize:4 /nologo /c /fast $*.F

.c.obj:
	cl /DCNS_ARCH_TYPE_WIN32 /c $*.c

#
# utilities
#

space:
	del *.obj

clean:
        del cns.exe
        del *.obj
#
