# Barnes-Hut Particle Simulation

This project simulates the gravitational interactions of particles using the Barnes-Hut algorithm. Below are the instructions to compile and run the simulation.

## Compilation

1. Create the build directory:

    ```sh
    mkdir build
    ```

2. Configure the project with CMake, enabling SDL and gprof profiling:

    ```sh
    cmake -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg -S . -B build -DUSE_SDL=ON
    ```

3. (bis) Configure the projec with CMake, enabling SDL and gprof as well as MPI and OpenMP:

    ```sh
    cmake -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg -S . -B build -DUSE_SDL=ON -DUSE_MPI=ON -DUSE_OPENMP=ON
    ```

4. Build the project:

    ```sh
    cmake --build build
    ```

## Running the Simulation

The binary for the simulation will be located in the `build` directory. You can run the simulation with the following command:

```sh
./build/particleSimulation <number_of_particles> <number_of_iterations> [sdl]
```

### Running the MPI+OpenMPhybrid version

The binaries created from the cmake compilation with the needed compilation flags need to be executed like the following:

```sh   
mpiexec -n <number_of_process> -x OMP_NUM_THREADS=<number_of_threads> ./build/particleSimulation <number_of_particles> <number_of_iterations>
```

You can also use the provided `start_mpi.sh` script that simplifies this command.

### Example

To run the simulation with 100 particles (in addition to the 9 default celestial bodies), 100000 iterations, and SDL graphical output:

```sh
./build/particleSimulation 100 100000 sdl
```

## Profiling

To generate a graphical gprof profiling output:

```sh
gprof ./build/particleSimulation | gprof2dot | dot -Tpng -o output.png
```

This will create a `output.png` file with the profiling information.
