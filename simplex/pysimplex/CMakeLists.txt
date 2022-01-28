include(GNUInstallDirs)

include(FetchContent)
FetchContent_Declare(
  pybind11
  URL https://github.com/pybind/pybind11/archive/refs/tags/v2.9.0.tar.gz
  URL_HASH SHA256=057fb68dafd972bc13afb855f3b0d8cf0fa1a78ef053e815d9af79be7ff567cb
)
FetchContent_MakeAvailable(pybind11)
pybind11_add_module(_simplex MODULE pysimplex.cpp)

target_link_libraries(_simplex PUBLIC simplex)

install(TARGETS _simplex DESTINATION .)

set(ignoreMe "${SKBUILD}")
