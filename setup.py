import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='asymmetree',
    version='0.1.0',
    author='David Schaller',
    author_email='david@0x002A.com',
    description='Simulation of phylogenetic trees, sequences and genomes.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/david-schaller/AsymmeTree',
    packages=setuptools.find_packages(exclude=['examples', 'examples.*',
                                               'resources', 'resources.*',
                                               'tests', 'tests.*',
                                               'validation', 'validation.*']),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    install_requires=[
        'numpy>=1.16.4',
        'scipy>=1.3.0',
        'matplotlib>=3.0',
        'networkx>=2.2',
   ],
)
