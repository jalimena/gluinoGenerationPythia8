# gluinoGenerationPythia8
files to create gluino r-hadrons for stopped particles experiment

```
tar xvfz pythia8306.tgz
cd pythia8306
./configure --with-root
```

maybe 
```source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.22.00/x86_64-centos7-gcc48-opt/bin/thisroot.sh```
or
```source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.22.00/x86_64-centos7-gcc48-opt/bin/thisroot.csh```
depending on your shell?

edit Makefile.inc to set:
```HEPMC2_USE=true```
and
```HEPMC3_USE=true```
So then it should look like my ```Makefile.inc``` included here

do:
```gmake```

copy gluinoGeneration.cc, HSCP_gluino_300_SLHA.spc, and Makefile to the examples/ directory

compile gluinoGeneration.cc with:
```make gluinoGeneration```

run gluinoGeneration.cc with:
```./gluinoGeneration```
you can do maybe:
```nohup ./gluinoGeneration >& out &```

it should just take a few seconds to generate 200 events