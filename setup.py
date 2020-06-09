import setuptools
import re

with open("README.md", "r") as fh:
    long_description = fh.read()

VERSIONFILE="CIAlign/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

setuptools.setup(
     name='cialign',
     version=verstr,
     author="Charlotte Tumescheit, Katy Brown",
     author_email="kab84@cam.ac.uk",
     description="Toolkit for cleaning and interpreting multiple sequence alignments",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/KatyBrown/CIAlign",
     packages=setuptools.find_packages(),
     install_requires=['matplotlib', 'numpy', 'ConfigArgParse', 'pillow'],
     scripts=['CIAlign/CIAlign'],
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent"]
 )
