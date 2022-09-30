# Lab5015SiPMModel
collection of programs to analyze SiPM pulse shape and compare them with the analytical model



### Login to your favourite machine
```sh
ssh username@hostname
```



### Fresh installation of the analysis package
```sh
export MYNAME=putYourNameHere  #use your name for development
mkdir $MYNAME
cd $MYNAME
git clone --recursive https://github.com/Lab5015/Lab5015SiPMModel
cd Lab5015SiPMModel
source scripts/setup.sh
make
make exe
```


### Package structure
This package is structured as follows
- Utilities are organized in the `interface` and `src` folders:
    - `interface`: contains the class or functions headers
    - `src`: contains the class or functions implementation
- `main`: contains the main analysis code that make use of the Utilities
- `cfg`: contains the config files which are used to pass parameters to the executables

After the compilation, each step of the analysis can be executed from the main folder `Lab5015Analysis` with a command like:
`./bin/executable cfg/configfile`



### Before running the analysis
The output filename and output location of each analysis step are defined in the cfg files. Before running the analysis make sure you have checked and if needed updated the relevant output paths in order not to overwrite the work of others.

Other than that, every time you login remember to source the setup script:
```
cd Lab5015SiPMModel
source scripts/setup.sh
```
