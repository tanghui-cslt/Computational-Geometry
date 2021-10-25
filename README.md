
## Get started with:

1. Downloading res
```bash
git clone --recursive https://github.com/tanghui-cslt/Computational-Geometry.git
```

* Windows

The core libigl functionality only depends on the C++ Standard Library and
Eigen.

To build all the examples in the tutorial, you can use the CMakeLists.txt in
the tutorial folder:
On my computer, 

```bash
cd tutorial
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -G "Visual Studio 15 2017 Win64" ../
```
and open `libigl_tutorials.sln` via visual stduio.
For different version of visual studio,  the instruction of cmake is different. such as, for visual studio 2015, you should use " Visual Studio 14 2015" to specify it.  So the last instruction should be modified by:

```bash
cd tutorial
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -G "Visual Studio 14 2015 Win64" ../
```

Notice: you must choose x64 compiler on windows.


* Linux

```bash
cd tutorial
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release  ../
make -j8
```

Notice: I just compile it on wsl2 (ubuntu 18.04) successfully, and I think it works on linux as well.

## Example

* 1-Conformal Mapping.

* Windows

you should set 1-Conformal_Mapping as Starup project firstly, and then run it.
![](./setting-1.png)

* Linux 

in `build` directory, 

```bash
./1-Conformal_Mapping
```

it's a conformal mapping, which maps a closed 2-dim surface to a ball conformally. 


* 2-Conformal Mapping.

* Windows

you should set 2-Harmonic-1-form as Starup project firstly, and then run it.


* Linux 

in `build` directory, 

```bash
./2-Harmonic-1-form
```

it's a conformal mapping, which maps a closed 2-dim surface to a ball conformally. 