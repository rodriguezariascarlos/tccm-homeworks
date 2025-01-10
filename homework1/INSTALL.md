# Installation Instructions

## Requirements
Before compiling and running the program, ensure the following dependencies are installed:

1. **TREXIO Library**: A library for handling quantum chemistry integrals.
2. **C Compiler**: Any standard C compiler such as GCC or Clang.

## Installation Steps

1. Install System Dependencies
   - On Ubuntu/Debian:
     ```bash
     sudo apt install libhdf5-dev
     ```
   - On macOS:
     ```bash
     brew install hdf5
     ```

2. Download and Install TREXIO Library
   - Download the latest release of TREXIO:
     [TREXIO 2.5.0](https://github.com/TREX-CoE/trexio/releases/download/v2.5.0/trexio-2.5.0.tar.gz)
   - Extract and install it:
     ```bash
     tar -xzf trexio-2.5.0.tar.gz
     cd trexio-2.5.0
     ./configure
     make
     sudo make install
     ```

3. Clone the Repository
   ```bash
   git clone <repository_url>
   cd <repository_name>
   ```

4. Compile the Program
   Use the following command to compile the source code:
   ```bash
   gcc -g -o main main.c -I/usr/local/include -L/usr/local/lib -ltrexio
   ```
   
5. Run the Program
   Execute the compiled binary by providing the path to a TREXIO file:
   ```bash
   ./main <molecule_file>
   ```

## Notes
- Ensure the TREXIO file provided as input contains valid molecular orbital data.
- For debugging or development, compile with `-g` to include debug information:
  ```bash
  gcc -g -o main main.c -I/usr/local/include -L/usr/local/lib -ltrexio 
  ```

