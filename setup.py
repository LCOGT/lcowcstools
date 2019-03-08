import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lcowcstools",
    version="0.0.1",
    author="Daniel Harbeck",
    author_email="dharbeck@lco.global",
    description="Toosl to model and fit SIP distortion coefficients in LCO and general FITS images.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/LCOGT/lcowcstools",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)