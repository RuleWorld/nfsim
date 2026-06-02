cmake_minimum_required(VERSION 3.5...3.27)

set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64" CACHE STRING "Build universal binary" FORCE)

project(NFsim)

set(CMAKE_BUILD_TYPE Release)

# Optional targets for embedding NFsim as a library or building the CLI.
option(NFSIM_BUILD_EXECUTABLE "Build the NFsim command line executable" ON)
option(NFSIM_BUILD_LIBRARY "Build NFsim as a reusable library target" OFF)

# Option to use ExprTk instead of muParser for expression parsing
option(NFSIM_USE_EXPRTK "Use ExprTk instead of muParser (requires C++17)" OFF)

if(NOT NFSIM_BUILD_EXECUTABLE AND NOT NFSIM_BUILD_LIBRARY)
    message(FATAL_ERROR "Enable at least one of NFSIM_BUILD_EXECUTABLE or NFSIM_BUILD_LIBRARY.")
endif()

set(SUB_DIRS
    src/nauty24
    src/NFutil
    src/NFutil/MTrand
    src/NFtest/transcription
    src/NFtest/tlbr
    src/NFtest/simple_system
    src/NFtest/scheduler
    src/NFtest/transformations
    src/NFtest/input
    src/NFtest/agentcell/cell
    src/NFtest/agentcell
    src/NFtest/util
    src/NFtest/tinyxml
    src/NFtest/nauty24
    src/NFtest/system
    src/NFtest/compartment
    src/NFtest/observable
    src/NFtest/mapping
    src/NFtest/molecule
    src/NFtest/moleculeType
    src/NFtest/complex
    src/NFtest/mappingSet
    src/NFscheduler
    src/NFreactions/transformations
    src/NFreactions/reactions
    src/NFreactions/reactantLists
    src/NFreactions/mappings
    src/NFreactions
    src/NFoutput
    src/NFinput
    src/NFinput/TinyXML
    src/NFinput/json
    src/NFfunction
    src/NFcore
    src/NFcore/reactionSelector
    src/NFcore/moleculeLists
)

# Conditionally add muParser or ExprTk to include path
if(NFSIM_USE_EXPRTK)
    list(APPEND SUB_DIRS src/NFfunction/exprtk)
else()
    list(APPEND SUB_DIRS src/NFfunction/muParser)
endif()

if(CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    set(BUILD_SHARED_LIBRARIES OFF)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR} include ${SUB_DIRS})

file(GLOB_RECURSE src_files  "${CMAKE_CURRENT_SOURCE_DIR}/src/*cpp")
file(GLOB_RECURSE c_src_files  "${CMAKE_CURRENT_SOURCE_DIR}/src/*c")

# Exclude muParser sources when using ExprTk
if(NFSIM_USE_EXPRTK)
    list(FILTER src_files EXCLUDE REGEX ".*/muParser/.*")
endif()

# Some branches contain both test_util.cpp and util_test.cpp defining the same
# NFtest_util::run() symbol. Prefer test_util.cpp and drop util_test.cpp to
# keep a single definition during link.
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/NFtest/util/util_test.cpp")
    list(FILTER src_files EXCLUDE REGEX "[/\\\\]src[/\\\\]NFtest[/\\\\]util[/\\\\]util_test\\.cpp$")
endif()

# The CLI entrypoint is built separately so the shared implementation can also
# be packaged as a library.
list(FILTER src_files EXCLUDE REGEX "[/\\\\]src[/\\\\]NFsim_main\\.cpp$")

set(SRC_FILES
    ${src_files}
    ${c_src_files}
)

set(NFSIM_OBJECT_TARGET "${PROJECT_NAME}_objects")
add_library(${NFSIM_OBJECT_TARGET} OBJECT ${SRC_FILES})

if(NFSIM_BUILD_LIBRARY)
    set_target_properties(${NFSIM_OBJECT_TARGET} PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

# Set C++ standard based on parser choice
if(NFSIM_USE_EXPRTK)
    set(NFSIM_CXX_STANDARD_FEATURE cxx_std_17)
    target_compile_definitions(${NFSIM_OBJECT_TARGET} PRIVATE NFSIM_USE_EXPRTK)
else()
    set(NFSIM_CXX_STANDARD_FEATURE cxx_std_11)
endif()

target_compile_features(${NFSIM_OBJECT_TARGET} PRIVATE ${NFSIM_CXX_STANDARD_FEATURE})

if(NOT MSVC)
    target_compile_options(${NFSIM_OBJECT_TARGET} PRIVATE -Wno-deprecated-declarations)
endif()

if(NFSIM_BUILD_LIBRARY)
    add_library(nfsim $<TARGET_OBJECTS:${NFSIM_OBJECT_TARGET}>)
    target_include_directories(nfsim PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_compile_features(nfsim PUBLIC ${NFSIM_CXX_STANDARD_FEATURE})

    if(NFSIM_USE_EXPRTK)
        target_include_directories(nfsim PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/NFfunction/exprtk")
        target_compile_definitions(nfsim PUBLIC NFSIM_USE_EXPRTK)
    endif()
endif()

if(NFSIM_BUILD_EXECUTABLE)
    add_executable(${PROJECT_NAME} src/NFsim_main.cpp $<TARGET_OBJECTS:${NFSIM_OBJECT_TARGET}>)
    target_compile_features(${PROJECT_NAME} PRIVATE ${NFSIM_CXX_STANDARD_FEATURE})

    if(NFSIM_USE_EXPRTK)
        target_compile_definitions(${PROJECT_NAME} PRIVATE NFSIM_USE_EXPRTK)
    endif()

    if(NOT MSVC)
        target_compile_options(${PROJECT_NAME} PRIVATE -Wno-deprecated-declarations)
    endif()
endif()
