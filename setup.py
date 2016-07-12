from setuptools import setup, find_packages

setup(
    name='mpmorph',
    version='0.1',
    packages=find_packages(),
    url='https://github.com/aykol/mpmorph',
    license='modified BSD',
    author='Muratahan Aykol',
    author_email='maykol@lbl.gov',
    description='',
    install_requires=['matmethods>=1.2.7', 'fireworks>=1.2.7', 'scipy>=0.17.1', 'matplotlib>=1.6']
)
