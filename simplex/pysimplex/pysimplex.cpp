#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <simplex.hpp>

#include <optional>
#include <ostream>
#include <sstream>

namespace py = pybind11;

PYBIND11_MODULE(_simplex, m) {
  py::class_<Simplex>(m, "Simplex",
    "Clifford circuit simulator")
    .def(py::init<unsigned>(),
        "Initialize a simulator with `n` qubits.",
        py::arg("n"))
    .def("__repr__",
        [](const Simplex& S) {
            std::stringstream ss;
            ss << "Simplex(" << std::endl << S << ")";
            return ss.str();
        })
    .def("X",
        &Simplex::X,
        "Apply an X gate to qubit `j`.",
        py::arg("j"))
    .def("Y",
        &Simplex::Y,
        "Apply a Y gate to qubit `j`.",
        py::arg("j"))
    .def("Z",
        &Simplex::Z,
        "Apply a Z gate to qubit `j`.",
        py::arg("j"))
    .def("H",
        &Simplex::H,
        "Apply an H gate to qubit `j`.",
        py::arg("j"))
    .def("S",
        &Simplex::S,
        "Apply an S gate to qubit `j`.",
        py::arg("j"))
    .def("Sdg",
        &Simplex::Sdg,
        "Apply an inverse S gate to qubit `j`.",
        py::arg("j"))
    .def("CX",
        &Simplex::CX,
        "Apply a CX gate to qubits `j` and `k`.",
        py::arg("j"), py::arg("k"))
    .def("CZ",
        &Simplex::CZ,
        "Apply a CZ gate to qubits `j` and `k`.",
        py::arg("j"), py::arg("k"))
    .def("MeasX",
        [](Simplex& S, unsigned j, std::optional<int> coin) {
            return S.MeasX(j, coin);
        },
        "Measure qubit `j` in the X basis."
        "\n\n"
        "If the outcome is non-deterministic, then the outcome is the value of"
        "`coin` if specified (must be 0 or 1), otherwise random.",
        py::arg("j"), py::arg("coin") = std::nullopt)
    .def("MeasY",
        [](Simplex& S, unsigned j, std::optional<int> coin) {
            return S.MeasY(j, coin);
        },
        "Measure qubit `j` in the Y basis."
        "\n\n"
        "If the outcome is non-deterministic, then the outcome is the value of"
        "`coin` if specified (must be 0 or 1), otherwise random.",
        py::arg("j"), py::arg("coin") = std::nullopt)
    .def("MeasZ",
        [](Simplex& S, unsigned j, std::optional<int> coin) {
            return S.MeasZ(j, coin);
        },
        "Measure qubit `j` in the Z basis."
        "\n\n"
        "If the outcome is non-deterministic, then the outcome is the value of"
        "`coin` if specified (must be 0 or 1), otherwise random.",
        py::arg("j"), py::arg("coin") = std::nullopt)
    .def("is_deterministic",
        &Simplex::is_deterministic,
        "Report whether all measurements are deterministic.");
}
