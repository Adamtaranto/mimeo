from setuptools import setup

pypi_classifiers = [
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    "Development Status :: 4 - Beta",
    "Environment :: Console",
    "Operating System :: OS Independent",
    'Intended Audience :: Science/Research',
    'Natural Language :: English',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    "Topic :: Software Development :: Libraries :: Python Modules",
    'License :: OSI Approved :: MIT License',
]

install_requires = [
    "pandas>=0.20.3",
    'biopython>=1.70',
]

desc = """Scan genomes for internally repeated sequences, elements which are repetitive in another species, or high-identity HGT candidate regions between species."""

setup(name='mimeo',
      version='0.1.0',
      description=desc,
      url='https://github.com/Adamtaranto/mimeo',
      author='Adam Taranto',
      author_email='adam.taranto@anu.edu.au',
      license='MIT',
      packages=['mimeo'],
      classifiers=pypi_classifiers,
      keywords=["Transposon","TE","WGA","LASTZ","Whole genome alignment","repeat","transposition"],
      install_requires=install_requires,
      include_package_data=True,
      zip_safe=False,
      entry_points={
        'console_scripts': [
            'mimeo-self=mimeo.run_self:main',
            'mimeo-x=mimeo.run_interspecies:main',
            'mimeo-map=mimeo.run_map:main',
        ],
    },
    )
