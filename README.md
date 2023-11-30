# LRS lib interface for matlab

This project includes a `cmake` project to compile the lrslib into a Matlab mex file.

## Dependencies

The only thing you will need is a `cmake` installation. This project includes its dependencies via `vcpkg`. To build it simply run

```
cmake -Bbuild -G Ninja
cmake --build build --target all
```

If you don't have `ninja` you can use your system's default generator as well, but in particular on Windows the output becomes significantly more intelligible using Ninja.