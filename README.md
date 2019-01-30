# CrystalPlan
Experiment planning tool for elastic neutron diffraction single-crystal experiments.

## Installing Locally for Development

The prerequisites for installing CrystalPlan are that you need the following
libraries available on your machine:
    
    - vtk
    - wxpython == 3.0
    - xtensor
    - xtensor-python
    - pybind11
    - Mayavi
    - gcc > 4.8 (c++14 support)

The easiest way to get all of these dependencies in one place is to use a conda
env. The code for this project can be installed using pip from the git repo:

``` 
pip install -e .  
```

This will install the required python dependencies and compile the C++ extensions.

## Running

CrystalPlan can be run by executing the following in the top level folder of
the git repo.

```
python crystalplan.py
```

### Running using Docker

The easiest way to run CrystalPlan is to grab the docker container and run it using xhost:

```
docker pull samueljackson92/crystal-plan
```

On OSX you can then run:

```
docker run --rm -it -e DISPLAY=docker.for.mac.localhost:0 --privileged -v /tmp/.X11-unix:/tmp/.X11-unix crystal-plan
```
The command for linux may differ a little in the arguments supplied to enable X11 forwarding from the container.

## Packaging

For packaging CrystalPlan for distribution to users and IDAaaS the PyInstaller
library is used. This tool finds all of the dependencies for the project and
packs them into a single folder that can be zipped and sent to the user.
PyInstaller supports all three operating systems and a different .spec file is
included in the repo for each one.

### CentOS / RHEL7
The easiest way is to grab a copy of the crystal-plan-package docker image.
Running this will generate a fully packaged version of CrystalPlan that can be
shipped to IDAaaS, isis-compute, or directly to end users.  

The docker container can be run as follows:

```
docker run -v $PWD:/CrystalPlan crystal-plan-dev
```

Where $PWD is the top level of your local crystal plan git repo. This will be
mounted to the container and then built for CentOS. The output of PyInstaller
will be directed into `./dist`. The packaged application can be run by
extracting the folder from the tar ball and executing the `crystalplan` binary.

### OSX

You'll need install all the relevant dependencies to a windows machine. The
easiest way I've found to do this is with conda where you can `conda install`
pretty much everything you need. You'll also need a compiler capable of
compiling at least C++14.

Then you're ready to `python setup.py build_ext` in the git repo to build the
C++ extensions. Once you've done this you should be able to run CrystalPlan on
OSX.

```
pythonw crystalplan.py
```

To package it for OSX you can run PyInstaller directly:

```
pyinstaller crystalplan-osx.spec -y
```

This will output the package application into the `dist` folder. To create a
.dmg for crystal plan you can use the `dmgbuild` tool (available on PyPI). Then
you can run:

```
dmgbuild -s ../settings.py "CrystalPlan" crystal-plan-osx.dmg
```

Which will output a dmg for CrystalPlan that can be shared.

### Windows

You'll need install all the relevant dependencies to a windows machine. The
easiest way I've found to do this is with conda where you can `conda install`
pretty much everything you need. You'll also need to install MS visual studio
>= 2015 and run the following commands in the command prompt prior to building
the C++ extensions for crystal plan.

```batch
"%VS140COMNTOOLS%\..\..\VC\vcvarsall.bat" x64
set DISTUTILS_USE_SDK=1
set MSSdk=1
```

Then you're ready to `python setup.py build_ext` in the git repo to build the
C++ extensions. Once you've done this you should be able to run CrystalPlan on
Windows:

```
python crystalplan.py
```

To package it for windows you can run PyInstaller directly:

```
pyinstaller crystalplan-windows.spec -y
```

The output can be found in the `dist` folder. The contents is self contained and
can be zipped and sent instrument scientists and other users.
