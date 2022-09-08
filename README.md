# RBFinterpolation
C++ code for RBF training and interpolation on a gridded Netcdf datasets.

The training is done automatically for each cell independantly.


## How to build

```
g++ RBF.cpp RBF.h read_input.cpp -I../include/ -I../include/armadillo_bits/  -llapack -lblas -lnetcdf -lgfortran
```

## Using

The code loolks for a `RBF_param.txt` in the current folder. 

### Training
```
#Header

centersfile = Centres_cascadia_3100.3_param.txt;

trainingfile = Samoa_mexico_Training_sub.nc?z;

RBFcoefffile=Samoa_cascadia_coeff.nc?3Dvar
#gammafile=Samoa_cascadia_coeff.nc?gamma
#inputfile = demoinput.txt
isdir=-1;
trainRBF = 1;
saveRBFcoeffs = 1;
interpRBF = 0;
#outputfile = SmallEQ.nc;
```

### Using

```
#Header

centersfile = Centres_cascadia_3100.3_param.txt;

#trainingfile = Samoa_mexico_Training_sub.nc?z;

RBFcoefffile=Samoa_cascadia_coeff.nc?3Dvar
gammafile=Samoa_cascadia_coeff.nc?gamma
inputfile = demoinput.txt
isdir=-1;
trainRBF = 0;
saveRBFcoeffs = 0;
interpRBF = 1;
outputfile = SmallEQ.nc;
```
