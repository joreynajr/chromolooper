from setuptools import setup
setup(
    name='chromolooper',
    version='0.1.0',
    author='Joaquin Reyna',
    author_email='joreyna@live.com',
    packages=['chromolooper'],
    scripts=[],
    url='',
    license='LICENSE.txt',
    description='Algebraic manipulation of chromosome loop data',
    long_description=open('README.md').read(),
    install_requires=["pybedtools"]
   )
