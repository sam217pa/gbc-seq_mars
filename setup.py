from setuptools import setup, find_packages

setup(name="phruscle",
      version="0.0.1",
      py_modules=["phruscle"],
      packages=find_packages('src'),
      package_dir={'': 'src'},
      install_requires=[
          'Click',
          'Biopython',
          'numpy',
          'pandas',
      ],
      entry_points='''
        [console_scripts]
        phruscle=phruscle:phruscle
    ''', )
