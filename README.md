# WellRadPy

WellRadPy ("Well Radius Python") is a Python package for calculating the radius
of influence and radius of investigation of pumping wells. These quantities
frequently need to be estimated in well hydraulics and aquifer testing (also
called well testing in petroleum engineering). In particular, calculating the
radius of investigation is critical for the design and the interpretation of
pumping tests. This package implements solutions for both drawdown and recovery
phases of pumping tests.

Note that radius of influence and radius of investigation are different
concepts. Furthermore, there are various ways by which each radius may be
precisely defined. This package implements most definitions.

## References

* Bresciani, E., Shandilya, R.N., Kang, P.K., Lee, S., 2020. Well radius of
influence and radius of investigation: What exactly are they and how to
estimate them? Journal of Hydrology, 583: 124646.
https://doi.org/10.1016/j.jhydrol.2020.124646.

* A paper focusing on radius of investigation during recovery is currently
under review.

## Usage

WellRadPy is simple and its use is straightforward. The ``/examples`` directory
in the [GitHub repository](https://github.com/etiennebresciani/wellradpy)
contains scripts where the functions for calculating radius of influence and
radius of investigation are exemplified.

Help can also be obtained by typing ``help(function_path)`` in the Python
console.

## Installation

### For simple use

1. From the command line, type ``pip install wellradpy``. This will download
the latest release of the package from the PyPI repository and install it.

### For development

1. Clone (or download) the sources from the
[GitHub repository](https://github.com/etiennebresciani/wellradpy).
2. From the project directory, type ``pip install -e .``.
This will install the current directory, so that you can both modify the
package and use it at the same time.

### Dependencies

WellRadPy depends on the Python packages NumPy and SciPy (>= 1.2.0). These
will be automatically installed with any of the above installation methods.

## Authors

* **Etienne Bresciani**
([etiennebresciani](https://github.com/etiennebresciani))

## License

This project is licensed under the MIT License &ndash; see the
[LICENSE.txt](LICENSE.txt) file for details.

## Acknowledgements

* Korea Ministry of Science and ICT through the Korea Research Fellowship
program of the National Research Foundation of Korea (2016H1D3A1908042)
* Korea Institute of Science and Technology (KIST) through the Future Research
Program (2E29660, 2E30510)
* Korea Ministry of Environment (MOE) through the Subsurface Environment
Management (SEM) Project of the Korea Environment Industry & Technology
Institute (KEITI) (2018002440006)
