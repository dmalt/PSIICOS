PSIICOS
=======

Toolbox for power and shift invariant imaging of coherent sources

Getting Started
---------------

### Prerequisites

#### OS
Installation script works only on LINUX. The only thing it does is
adding PSIICOS folder location to the MATLAB path. It can be easily done manually.

The rest is cross-platform

#### Packages

PSIICOS depends on:
 * [utils_psiicos](https://github.com/dmalt/utils_psiicos.git)

### Installing

In LINUX run in terminal:

```bash
git clone https://github.com/dmalt/PSIICOS.git && cd $_
./install.sh
```

In WINDOWS get the package by downloading .zip archive from github,
extract the files and add path to package folder from inside MATLAB.

To test the installation open matlab, change directory to 
`PSIICOS/tests` and run

```matlab
check_installation
```

It should print

```matlab
INSTALLATION SUCCESFUL
```

If it didn't, check paths and script install.sh
