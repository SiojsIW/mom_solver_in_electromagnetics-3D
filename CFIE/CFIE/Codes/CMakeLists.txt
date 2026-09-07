cmake_minimum_required(VERSION 3.20)
project(fortran_practice Fortran)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# 扫描所有 .f90
file(GLOB ALL_F90 "*.f90")

# 分类：纯模块（无PROGRAM） vs 主程序文件（有PROGRAM）
set(MODULE_ONLY "")
set(MAIN_FILES "")

foreach(file ${ALL_F90})
    file(STRINGS ${file} has_program REGEX "^[ \t]*PROGRAM[ \t]+")
    if(has_program)
        list(APPEND MAIN_FILES ${file})
    else()
        list(APPEND MODULE_ONLY ${file})
    endif()
endforeach()

# 纯模块编译成公共库，所有主程序自动链接
if(MODULE_ONLY)
    add_library(common_modules STATIC ${MODULE_ONLY})
endif()

# 每个主程序单独生成可执行文件
foreach(file ${MAIN_FILES})
    get_filename_component(target_name ${file} NAME_WE)
    string(REPLACE "." "_" target_name ${target_name})
    add_executable(${target_name} ${file})
    if(MODULE_ONLY)
        target_link_libraries(${target_name} common_modules)
    endif()
endforeach()