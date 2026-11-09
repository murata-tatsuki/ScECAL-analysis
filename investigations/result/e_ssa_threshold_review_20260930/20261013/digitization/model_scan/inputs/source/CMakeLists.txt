cmake_minimum_required(VERSION 3.0)

# 设置项目名称
project(MCAnalysis)

# 设置C++标准
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 包含ROOT的头文件目录
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)

# 包含头文件目录
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})

# 查找所有源文件
file(GLOB SOURCES "src/*.cxx")

# 添加可执行文件
add_executable(MCAnalysis main.cxx ${SOURCES})

# 链接ROOT库
target_link_libraries(MCAnalysis ${ROOT_LIBRARIES})
