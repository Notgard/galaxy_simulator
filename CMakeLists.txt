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
option(USE_MPI "Enable MPI for parallel execution" ON)  # Still keep this option
option(USE_OPENMP "Enable OpenMP for parallel execution" ON)

set(CUSTOM_COMPILE_FLAGS "-pg -g -gdwarf-3 -O0 -march=native")
string(REPLACE " " ";" CUSTOM_COMPILE_FLAGS_LIST ${CUSTOM_COMPILE_FLAGS})

# Source files (conditionally include MPI versions)
set(SOURCES
    src/main.cpp      # default main file
    src/Simulation.cpp # default Simulation file
)

# If MPI is enabled, use the MPI-specific files
if(USE_MPI)
    list(TRANSFORM SOURCES REPLACE src/main.cpp src/mpi_main.cpp)  # replace main.cpp with mpi_main.cpp
    list(APPEND SOURCES src/SimulationMPI.cpp)
    #list(TRANSFORM SOURCES REPLACE src/Simulation.cpp src/SimulationMPI.cpp)  # replace Simulation.cpp with SimulationMPI.cpp
endif()

# If SDL is enabled, use the SDL-specific files
if(USE_SDL)
    list(TRANSFORM SOURCES REPLACE src/main.cpp src/sdl_main.cpp)  # replace main.cpp with mpi_main.cpp
    list(APPEND SOURCES src/SimulationSDL.cpp)
endif()

set(HEADERS
    includes/SDLUtils.h
    includes/Box.h
    includes/Config.h
    includes/Quadtree.h
    includes/Simulation.h
    includes/SimulationMPI.h
    includes/SimulationSDL.h
    includes/Vector2d.h
    includes/Particle.h 
)

include_directories(includes)

# Add executable
add_executable(particleSimulation ${SOURCES} ${HEADERS})

target_include_directories(${CMAKE_PROJECT_NAME}
    PRIVATE ${CMAKE_SOURCE_DIR}/includes
)

target_compile_options(${CMAKE_PROJECT_NAME} PRIVATE ${CUSTOM_COMPILE_FLAGS_LIST})

# Find SDL2
find_package(SDL2 REQUIRED)
message("SDL version: " ${SDL2_VERSION})
include_directories(${SDL2_INCLUDE_DIRS})

if(USE_SDL)
    target_compile_definitions(${CMAKE_PROJECT_NAME} PRIVATE USE_SDL)

    # Find SDL_ttf
    find_package(SDL2_ttf REQUIRED)
    message("SDL_ttf found")
    target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE SDL2::Main SDL2::TTF)
endif()

# Find GLM
find_package(glm REQUIRED)

target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE glm::glm)

# MPI Support
if(USE_MPI)
    find_package(MPI REQUIRED)
    if(MPI_FOUND)
        message("MPI found: " ${MPI_CXX_COMPILER})
        target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${MPI_INCLUDE_PATH})
        target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE MPI::MPI_CXX)
        target_compile_definitions(${CMAKE_PROJECT_NAME} PRIVATE USE_MPI)
    else()
        message(FATAL_ERROR "MPI is required but not found.")
    endif()
endif()

# OpenMP Support
if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_CXX_FOUND)
        message("OpenMP found")
        target_compile_options(${CMAKE_PROJECT_NAME} PRIVATE -fopenmp)
        target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
        target_compile_definitions(${CMAKE_PROJECT_NAME} PRIVATE USE_OPENMP)
    endif()
endif()

# Add source files
target_sources(${CMAKE_PROJECT_NAME}
    PRIVATE ${sources}
)
