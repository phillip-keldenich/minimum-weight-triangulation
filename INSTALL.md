This code uses a combination of `CMake` and `conan` to build.
With one exception, all dependencies are available on official
`conan` repositories; however, that exception (`gurobi`) needs to be installed
before attempting to build this code.
For that purpose, go to the `conan_pkg/gurobi_public` directory and run

```
conan create . --user ibralg --channel develop -s build_type=Release --build=missing
```


Afterwards, go back to the directory containing this `INSTALL.md` and
run
```
conan install . -s build_type=Release --build=missing
```

to install the dependencies and setup the build environment.
Afterwards, you should be able to build the code using CMake.
Then, first run 
```
cmake -S . -B build/Release -DCMAKE_BUILD_TYPE=Release
```
followed by 
```
cmake --build build/Release 
```
to build the software (which will end up below `build/Release`).
