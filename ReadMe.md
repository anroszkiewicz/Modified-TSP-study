## Modified TSP study

This repository contains a comparison of multiple methods that can be used to solve a modified TSP problem, in which we construct two separate paths, each containg half of the vertices.

The test instances are taken from [TSPLib](http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/) and the results are visualized using [matplotlib-cpp](https://github.com/lava/matplotlib-cpp/blob/master/matplotlibcpp.h).

### Authors

Mikołaj Nowak, Anna Roszkiewicz

Intelligent optimization methods course

### Running the code

To run the code with `matplotlib-cpp`, you need to have `python3-dev` installed. Then, include links to Python's headers and libraries when compiling.

Windows:

```bash
g++ -std=c++11 main.cpp -o myplot -DWITHOUT_NUMPY -I"C:\Users\user\AppData\Local\Programs\Python\Python311\Include" -L"C:\Users\user\AppData\Local\Programs\Python\Python311\libs" -lpython311
```

Linux:

```bash
g++ -std=c++11 main.cpp -o myplot -DWITHOUT_NUMPY -I/usr/include/python3.12 -lpython3.12
```
