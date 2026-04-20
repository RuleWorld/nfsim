cmake_minimum_required(VERSION 3.5...3.27)

set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64" CACHE STRING "Build universal binary" FORCE)

project(NFsim)

set(CMAKE_BUILD_TYPE Release)

set(SUB_DIRS
    src/nauty24
    src/NFutil
    src/NFutil/MTrand
    src/NFtest/transcription
    src/NFtest/tlbr
    src/NFtest/simple_system
    src/NFtest/scheduler
    src/NFtest/transformations
    src/NFtest/agentcell/cell
    src/NFtest/agentcell
    src/NFtest/util
    src/NFtest/system
    src/NFtest/mapping
    src/NFtest/molecule
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
    src/NFfunction/muParser
    src/NFfunction
    src/NFcore
    src/NFcore/reactionSelector
    src/NFcore/moleculeLists
)

if(CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    set(BUILD_SHARED_LIBRARIES OFF)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR} include ${SUB_DIRS})

file(GLOB_RECURSE src_files  "${CMAKE_CURRENT_SOURCE_DIR}/src/*cpp")
file(GLOB_RECURSE c_src_files  "${CMAKE_CURRENT_SOURCE_DIR}/src/*c")

# Some branches contain both test_util.cpp and util_test.cpp defining the same
# NFtest_util::run() symbol. Prefer test_util.cpp and drop util_test.cpp to
# keep a single definition during link.
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/NFtest/util/util_test.cpp")
    list(FILTER src_files EXCLUDE REGEX "[/\\\\]src[/\\\\]NFtest[/\\\\]util[/\\\\]util_test\\.cpp$")
endif()

set(SRC_FILES
    ${src_files}
    ${c_src_files}
)

add_executable(${PROJECT_NAME} ${SRC_FILES})

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_11)

if(NOT MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE -Wno-deprecated-declarations)
endif()
