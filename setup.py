import os
from setuptools import setup, find_packages
from src.version import __version__

cur_path = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(cur_path, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
      name='Swave.py',
      version=__version__,

      description='SV/CSV callers on assemblies',
      long_description=long_description,

      url='https://github.com/songbowang125/',

      author='Songbo Wang',
      author_email='songbowang125@163.com',

      license='GPLv3',
      classifiers=[
      'Operating System :: POSIX :: Linux',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
      'Programming Language :: Python :: 3.7'
      ],

      keywords=['Swave', 'Structural variants(SV)', 'Complex structural variants', 'Assembly', 'Pangenome'],

      packages = ['src', 'src/pack_dotplot', 'src/pack_graph', 'src/pack_model', 'src/pack_sv'],
      data_files = [("", ["LICENSE"])],

      zip_safe=False,

      scripts=['Swave.py'],
      )