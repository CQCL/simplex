from skbuild import setup

setup(
    name="pysimplex",
    version="0.0.1",
    description="A fast Clifford circuit simulator",
    long_description="Based on https://arxiv.org/abs/2109.08629. Very alpha! Known bugs.",
    author="Alec Edgington",
    author_email="alec.edgington@cambridgequantum.com",
    license="Apache 2",
    url="https://github.com/CQCL/simplex",
    python_requires=">=3.8",
    classifiers=[
        "Environment :: Console",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: POSIX :: Linux",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
    ],
    packages=["pysimplex"],
    cmake_install_dir="pysimplex",
)
