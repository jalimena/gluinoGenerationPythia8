# gluinoGenerationPythia8
files to create gluino r-hadrons for stopped particles experiment

```
tar xvfz pythia8306.tgz
cd pythia8306
./configure --with-root
```

do
```source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.22.00/x86_64-centos7-gcc48-opt/bin/thisroot.sh```

or
```source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.22.00/x86_64-centos7-gcc48-opt/bin/thisroot.csh```
depending on your shell

edit Makefile.inc to set:
```HEPMC2_USE=true```
and
```HEPMC3_USE=true```
So then it should look like my ```Makefile.inc``` included here

do:
```gmake```

copy gluinoGeneration.cc, data/HSCP_gluino_300_SLHA.spc, and Makefile to the examples/ directory

compile gluinoGeneration.cc with:
```make gluinoGeneration```

run gluinoGeneration.cc with:
```./gluinoGeneration```
1st argument is gluino mass [GeV], 2nd is the width of the detector [m], and 3rd is detector position (0 for above collision point and 1 for diagonal to collision point). e.g. for 5GeV gluino, 2x2m detector above the collision point:
```nohup ./gluinoGeneration 5 2 0 >& out &```

it should just take a few seconds to generate 200 events

To run plotHist.cc open ROOT. do:
```.L plotHist.cc```
```.x plotHist.cc("2","0")```
where the 1st argument is the width of the detector [m], and the 2nd is detector postion