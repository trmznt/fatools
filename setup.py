import os

from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.md')) as f:
    README = f.read()
with open(os.path.join(here, 'CHANGES.md')) as f:
    CHANGES = f.read()

requires = [
    'numpy',
    'scipy',
    'matplotlib',
    'pyyaml',
    'pandas',
    ]

setup(name='fatools',
      version='0.7',
      description='fatools',
      long_description=README + '\n\n' + CHANGES,
      classifiers=[
        "Programming Language :: Python",
        ],
      author='Hidayat Trimarsanto',
      author_email='anto@eijkman.go.id',
      url='',
      keywords='dna fragment-analysis',
      packages=find_packages(),
      include_package_data=True,
      zip_safe=False,
      install_requires=requires,
      tests_require=requires,
      test_suite="fatools",
      entry_points="""\
      [console_scripts]
      fatools = fatools.scripts.run:main
      """,
      )
