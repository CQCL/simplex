include(GNUInstallDirs)

include(FetchContent)
FetchContent_Declare(
  pybind11
  URL https://github.com/pybind/pybind11/archive/refs/tags/v2.11.1.tar.gz
  URL_HASH SHA256=d475978da0cdc2d43b73f30910786759d593a9d8ee05b1b6846d1eb16c6d2e0c
)
FetchContent_MakeAvailable(pybind11)
pybind11_add_module(_simplex MODULE pysimplex.cpp)

target_link_libraries(_simplex PUBLIC simplex)

install(TARGETS _simplex DESTINATION .)

set(ignoreMe "${SKBUILD}")
