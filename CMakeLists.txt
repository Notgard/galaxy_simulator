cmake_minimum_required(VERSION 3.10)
project("particleSimulation"
    DESCRIPTION "N-Body simulation with SDL2 visualization"
    LANGUAGES C CXX
)

message("Project name: " ${CMAKE_PROJECT_NAME})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Options
option(USE_SDL "Enable SDL for visualization" ON)
option(USE_MPI "Enable MPI for parallel execution" ON)
option(USE_OPENMP "Enable OpenMP for parallel execution" ON)

# Custom flags
set(DEBUG_COMPILE_FLAGS "-O0")
set(OPTIMIZATION_COMPILE_FLAGS "-Ofast -march=znver4 -m64 -ffast-math -ffp-contract=fast -mavx512f")

# Apply optimization flags if not in debug mode
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    string(REPLACE " " ";" CUSTOM_COMPILE_FLAGS_LIST ${DEBUG_COMPILE_FLAGS})
    message("Compiler flags: ${DEBUG_COMPILE_FLAGS}")
else()
    string(REPLACE " " ";" CUSTOM_COMPILE_FLAGS_LIST ${OPTIMIZATION_COMPILE_FLAGS})
    message("Compiler flags: ${OPTIMIZATION_COMPILE_FLAGS}")
endif()

set(COMMON_SOURCES
    src/Simulation.cpp
)

set(COMMON_HEADERS
    includes/SDLUtils.h
    includes/Box.h
    includes/Config.h
    includes/Quadtree.h
    includes/Simulation.h
    includes/SimulationMPI.h
    includes/SimulationSDL.h
    includes/SimulationMPISDL.h
    includes/Vector2d.h
    includes/Particle.h
)

include_directories(includes)

# --- Non-MPI Version ---
if(NOT USE_MPI)
    set(NON_MPI_SOURCES ${COMMON_SOURCES} src/main.cpp)
    
    if(USE_SDL)
        list(APPEND NON_MPI_SOURCES src/SimulationSDL.cpp src/sdl_main.cpp)
    endif()

    add_executable(particleSimulation ${NON_MPI_SOURCES} ${COMMON_HEADERS})

    target_include_directories(particleSimulation PRIVATE ${CMAKE_SOURCE_DIR}/includes)
    target_compile_options(particleSimulation PRIVATE ${CUSTOM_COMPILE_FLAGS_LIST})

    # SDL Support
    if(USE_SDL)
        find_package(SDL2 REQUIRED)
        if(SDL2_FOUND)
            message("SDL version: " ${SDL2_VERSION})
            target_include_directories(particleSimulation PRIVATE ${SDL2_INCLUDE_DIRS})
            target_compile_definitions(particleSimulation PRIVATE USE_SDL)

            find_package(SDL2_ttf REQUIRED)
            message("SDL_ttf found")
            target_link_libraries(particleSimulation PRIVATE SDL2::Main SDL2::TTF)
        else()
            message(FATAL_ERROR "SDL2 not found")
        endif()
    endif()
endif()

# --- MPI Version ---
if(USE_MPI)
    set(MPI_SOURCES ${COMMON_SOURCES} src/mpi_main.cpp src/SimulationMPI.cpp)

    if(USE_SDL)
        list(APPEND MPI_SOURCES src/SimulationMPISDL.cpp src/SimulationSDL.cpp)
        list(TRANSFORM MPI_SOURCES REPLACE "src/mpi_main.cpp" "src/mpi_sdl_main.cpp")
    endif()

    add_executable(particleSimulationMPI ${MPI_SOURCES} ${COMMON_HEADERS})

    target_include_directories(particleSimulationMPI PRIVATE ${CMAKE_SOURCE_DIR}/includes)
    target_compile_options(particleSimulationMPI PRIVATE ${CUSTOM_COMPILE_FLAGS_LIST})

    find_package(MPI REQUIRED)
    if(MPI_FOUND)
        message("MPI found: " ${MPI_CXX_COMPILER})
        target_include_directories(particleSimulationMPI PRIVATE ${MPI_INCLUDE_PATH})
        target_link_libraries(particleSimulationMPI PRIVATE MPI::MPI_CXX)
        target_compile_definitions(particleSimulationMPI PRIVATE USE_MPI)
    else()
        message(FATAL_ERROR "MPI is required but not found.")
    endif()

    # SDL Support in MPI Version
    if(USE_SDL)
        find_package(SDL2 REQUIRED)
        if(SDL2_FOUND)
            message("SDL version: " ${SDL2_VERSION})
            target_include_directories(particleSimulationMPI PRIVATE ${SDL2_INCLUDE_DIRS})
            target_compile_definitions(particleSimulationMPI PRIVATE USE_SDL)

            find_package(SDL2_ttf REQUIRED)
            message("SDL_ttf found")
            target_link_libraries(particleSimulationMPI PRIVATE SDL2::Main SDL2::TTF)
        else()
            message(FATAL_ERROR "SDL2 not found")
        endif()
    endif()
endif()

# OpenMP Support
if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_CXX_FOUND)
        message("OpenMP found")
        
        if(TARGET particleSimulation)
            target_compile_options(particleSimulation PRIVATE -fopenmp)
            target_link_libraries(particleSimulation PRIVATE OpenMP::OpenMP_CXX)
            target_compile_definitions(particleSimulation PRIVATE USE_OPENMP)
        endif()

        if(TARGET particleSimulationMPI)
            target_compile_options(particleSimulationMPI PRIVATE -fopenmp)
            target_link_libraries(particleSimulationMPI PRIVATE OpenMP::OpenMP_CXX)
            target_compile_definitions(particleSimulationMPI PRIVATE USE_OPENMP)
        endif()
    endif()
endif()
