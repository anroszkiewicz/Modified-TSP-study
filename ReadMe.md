# README: Compiling C++ with Python Integration Using `g++`

This guide explains how to compile a C++ program that integrates with Python using `g++`. The example provided shows how to link your C++ code with Python libraries like `matplotlib` for visualization.

### Command Breakdown:

```bash
g++ -std=c++11 main.cpp -o myplot -DWITHOUT_NUMPY -I"C:\Users\user\AppData\Local\Programs\Python\Python311\Include" -L"C:\Users\user\AppData\Local\Programs\Python\Python311\libs" -lpython311
```

- **`g++`**: The GNU C++ compiler.
- **`-std=c++11`**: Specifies C++11 standard.
- **`main.cpp`**: Source file to compile.
- **`-o myplot`**: Output executable name (`myplot`).
- **`-DWITHOUT_NUMPY`**: Preprocessor macro (can exclude numpy code).
- **`-I`**: Include Python's headers (adjust path as needed).
- **`-L`**: Link to Python libraries (adjust path as needed).
- **`-lpython311`**: Links the Python 3.11 library.
