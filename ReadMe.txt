g++ -I/Users/nushamehmanesh/SSDR-2012/ svg_parser.cpp 
g++ -I/usr/local/Cellar/eigen/3.3.4/include/eigen3 -I/Users/nushamehmanesh/SSDR-2012/ svg_parser.cpp
g++ -I/usr/local/Cellar/eigen/3.3.4/include/eigen3 -I/Users/nushamehmanesh/SSDR-2012/ -o svg_parser svg_parser.cpp
./svg_parser ./svg_files/walk1-v2.svg
g++ -std=c++11 -I/usr/local/Cellar/eigen/3.3.4/include/eigen3 -I/Users/nushamehmanesh/SSDR-2012/ -o svg_parser svg_parser.cpp
./svg_parser 3 ./svg_files/walk1-v2.svg 2 ./svg_files/walk1-v2.svg ./svg_files/walk1-v2.svg


g++ -std=c++11 -I/usr/local/Cellar/eigen/3.3.4/include/eigen3 -I/Users/nushamehmanesh/SSDR-2012/ -I/usr/local/Cellar/gsl -o svg_parser svg_parser.cpp ssdr.cpp

./svg_parser 4 ./svg_files/walk1-v2.svg 2 ./svg_files/walk1-v2.svg ./svg_files/walk1-v2.svg

