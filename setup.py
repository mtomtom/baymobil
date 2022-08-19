from setuptools import setup, find_packages

VERSION = '0.1.5' 
DESCRIPTION = 'Bayesian analysis of heterograft data'
LONG_DESCRIPTION = 'A function to evaluate heterograft data based on error rate distributions taken from homograft data'

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="baymobil", 
        version=VERSION,
        author="Melissa Tomkins",
        author_email="melissa.tomkins@jic.ac.uk",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=["numpy","scipy","pandas"], # add any additional packages that 
        # needs to be installed along with your package.
        
        keywords=['python', 'first package'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
        ]
)