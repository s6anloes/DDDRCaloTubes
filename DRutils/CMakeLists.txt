# project_root/utils/CMakeLists.txt
set(CMAKE_CXX_STANDARD 17)

# Define the library
add_library(DRutils DRutils.cpp)

# Specify include directories for the library
target_include_directories(DRutils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(DRutils DD4hep::DDCore ROOT::Core )
# Set position-independent code property
set_target_properties(DRutils PROPERTIES POSITION_INDEPENDENT_CODE ON)