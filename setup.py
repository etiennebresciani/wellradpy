import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="wellradpy",
    version="1.0.0",
    author="Etienne Bresciani",
    description="A small Python package to calculate the radius of influence and radius of investigation of wells.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/etiennebresciani/wellradpy",
    keywords="groundwater wells hydraulics",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "numpy",
        "scipy>=1.2.0"],
)
