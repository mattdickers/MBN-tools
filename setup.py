from setuptools import setup, find_packages

setup(
    name="MBN-tools",
    version="1.0.0",
    author="Matthew Dickers",
    author_email="mattdickers@gmail.com",
    description="Useful tools for use with MBN Explorer simulations, files, and their analysis",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/mattdickers/MBN-tools",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
          'numpy',
          'mdtraj',
          'scipy'
      ],
)
