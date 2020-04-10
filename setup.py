import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='cialign',  
     version='1.0.4',
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
