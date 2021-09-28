import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='asymmetree',
    version='2.1.0',
    author='David Schaller',
    author_email='sdavid@bioinf.uni-leipzig.de',
    description='Simulation of phylogenetic trees, sequences and genomes.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/david-schaller/AsymmeTree',
    packages=setuptools.find_packages('src'),
    package_dir={'': 'src'},
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
    install_requires=[
        'numpy>=1.16.4',
        'scipy>=1.3.0',
        'matplotlib>=3.0',
        'networkx>=2.2',
        'tralda>=1.0.0',
   ],
)
