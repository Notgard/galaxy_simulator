cmake_minimum_required(VERSION 3.10)
project("particleSimulation"
    DESCRIPTION "N-Body simulation with SDL2 visualisation"
    LANGUAGES C CXX
)

message("Project name: " ${CMAKE_PROJECT_NAME})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Option to enable/disable SDL
option(USE_SDL "Enable SDL for visualization" ON)

set(CUSTOM_COMPILE_FLAGS "-pg -g -gdwarf-3 -O0 -march=native")  # Example: Optimization + warnings
string(REPLACE " " ";" CUSTOM_COMPILE_FLAGS_LIST ${CUSTOM_COMPILE_FLAGS})

# Add source files
set(SOURCES
    src/main.cpp
    src/Simulation.cpp
)

set(HEADERS
    includes/Box.h
    includes/Config.h
    includes/Quadtree.h
    includes/Simulation.h
    includes/Vector2d.h
    includes/Particle.h
)

# Add include directories
include_directories(includes)

# Add executable
add_executable(particleSimulation ${SOURCES} ${HEADERS})

target_include_directories(${CMAKE_PROJECT_NAME}
    PRIVATE
        ${CMAKE_SOURCE_DIR}/includes
)

target_compile_options(${CMAKE_PROJECT_NAME} PRIVATE ${CUSTOM_COMPILE_FLAGS_LIST})

# Find SDL2
find_package(SDL2 REQUIRED)
message("SDL version: " ${SDL2_VERSION})
include_directories(${SDL2_INCLUDE_DIRS})

if(USE_SDL)
    target_compile_definitions(${CMAKE_PROJECT_NAME} PRIVATE USE_SDL)
endif()

#Find GLM
find_package(glm REQUIRED)

# SDL2::SDL2main may or may not be available. It is e.g. required by Windows GUI applications
if(TARGET SDL2::SDL2main)
    # It has an implicit dependency on SDL2 functions, so it MUST be added before SDL2::SDL2 (or SDL2::SDL2-static)
    target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE SDL2::SDL2main)
endif()

# Link to the actual SDL2 library. SDL2::SDL2 is the shared SDL library, SDL2::SDL2-static is the static SDL libarary.
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE SDL2::SDL2 glm::glm)

# Add source files to the target
target_sources(${CMAKE_PROJECT_NAME}
    PRIVATE ${sources}
)
