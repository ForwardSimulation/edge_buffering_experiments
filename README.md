# Edge buffering for tree sequence recording.

Working out kinks in the idea... :)

# Python dependencies

* tskit
* attrs
* seaborn
* pandas

Versions used here for the important things:

* Python 3.8
* tskit 2.3 (pip3 installed from PyPi)
* attrs 19.3.0 (pip3 installed from PyPi)

# C++ dependencies

* A C++ compiler.
* boost's program options library, ``sudo apt install libboost-program-options-dev``
* GSL, `sudo apt install libgsl-dev`
* cmake, `sudo apt install cmake`

The tool chain used here is:

* Pop! OS 20.04
* g++ 9.3
* boost 1.71
* gsl 2.5
* cmake 3.16.3

All packages are the OS versions.

## Building the C++ example

```sh
cmake -Bbuild -DCMAKE_BUILD_TYPE=Release
cd build
make VERBOSE=1
```

Substitute `Debug` for `Release` to disable optimizations and add `-g`.
