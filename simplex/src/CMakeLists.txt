add_library(simplex
    simplex.cpp
    A_matrix.cpp
    Q_matrix.cpp)

target_include_directories(simplex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
