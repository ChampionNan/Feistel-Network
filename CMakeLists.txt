cmake_minimum_required(VERSION 3.6)
# 定義工程名稱
project(Enc)

set(SOURCE_FILES enc.cpp)
# Header file path
set(INC_DIR /usr/local/include)
# dynamic link library
set(LINK_DIR /usr/local/lib)
# 定義依賴動態連結庫

include_directories(${INC_DIR})

link_directories(${LINK_DIR})

link_libraries(mbedtls mbedcrypto mbedx509)

add_executable(enc ${SOURCE_FILES})

target_link_libraries(enc mbedtls)
