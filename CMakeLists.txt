cmake_minimum_required(VERSION 3.24)

project(RTCamp8 CXX)

set(TBB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/oneTBB" CACHE PATH "" FORCE)
set(TBB_INCLUDE_DIR "${TBB_ROOT}/include" CACHE PATH "" FORCE)
set(TBB_SANITIZE "" CACHE STRING "" FORCE)
set(TBB_TEST OFF CACHE BOOL "" FORCE)
add_subdirectory(ThirdParty/oneTBB)

set(OIDN_APPS OFF CACHE BOOL "" FORCE)
set(OIDN_STATIC_LIB ON CACHE BOOL "" FORCE)
add_subdirectory(ThirdParty/OpenImageDenoise)

set(XITILS_BUILD_SANDBOX OFF CACHE BOOL "" FORCE)
add_subdirectory(Xitils)

set(CMAKE_CXX_FLAGS_RELEASE "/MT")
set(CMAKE_CXX_FLAGS_DEBUG "/MTd")
enable_language(CXX)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Zi -DDEBUG")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON) 
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
set(CMAKE_CXX_FLAGS "/openmp")
add_definitions(-D_UNICODE)
add_definitions(-DNOMINMAX)
add_definitions(/MP)
add_definitions(/EHsc)

set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Data")

add_definitions(/EHsc)
add_executable(RTCamp8 main.cpp)
target_link_libraries(RTCamp8 PRIVATE
	Xitils
	OpenImageDenoise
	)
target_include_directories(RTCamp8 PRIVATE
	${XITILS_INCLUDE_DIR}
	)
set_target_properties(RTCamp8 PROPERTIES
   VS_DEBUGGER_WORKING_DIRECTORY ${DATA_DIR}
)

