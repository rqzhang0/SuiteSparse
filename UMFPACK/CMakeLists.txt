cmake_minimum_required(VERSION 3.0)
project(umfpack)

# 设置源代码目录
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Source)

# 查找所有的 .cpp 文件
file(GLOB CPP_SOURCES ${SOURCE_DIR}/*.cpp)

# 查找所有的 .c 文件
file(GLOB C_SOURCES ${SOURCE_DIR}/*.c)

# 添加动态库目标
add_library(umfpack SHARED ${CPP_SOURCES} ${C_SOURCES})

# 设置编译选项
target_compile_options(umfpack PRIVATE
    -Wall
    -Wextra 
		-fexceptions
		-fPIC
		-fopenmp
		-DDINT
    # 其他编译选项...
)

# 设置链接选项
target_link_libraries(umfpack PRIVATE
amd suitesparseconfig ccolamd camd metis rt cholmod m
)

target_link_directories(umfpack PRIVATE /home/rh/CODE/SuiteSparse/lib /usr/lib/x86_64-linux-gnu/blas /usr/lib/x86_64-linux-gnu/lapack)

target_include_directories(umfpack PRIVATE /home/rh/CODE/SuiteSparse/include /home/rh/CODE/SuiteSparse/AMD/Include /home/rh/CODE/SuiteSparse/SuiteSparse_config)


