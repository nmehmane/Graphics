All created svg files are located in the svg_files directory.

This Program uses the SVG++ header only library which is already included for your convinience.
SVG++ makes use of the boost library.

## Build instructions:

In the project directory:

```
mkdir build
cd build
cmake ..
make -j8
cd bin
```


## Running

You can run the program with:

`./svg_parser <number of handles> <svg rest-pose> <number of frames> <svg pose 1> <svg pose2> ...`

Example:

` ./svg_parser 4 ./svg_files/walk1-v2.svg 1 ../../svg_files/walk1-v2.svg`

`./svg_parser 2 ./svg_files/sample2-input-1.svg 1 ./svg_files/sample2-input-2.svg`


## Misc

Building without cmake:

```
//g++ -std=c++11 -I/usr/local/Cellar/eigen/3.3.4/include/eigen3 -I/Users/nushamehmanesh/SSDR-2012/ -I/usr/local/Cellar/gsl -o svg_parser svg_parser.cpp ssdr.cpp
```
