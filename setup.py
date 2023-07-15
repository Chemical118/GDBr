from gdbr.version import get_version
from setuptools import setup


with open('README.md', 'r') as f:
    long_description = f.read()

setup(name='GDBr',
      version=get_version(),
      author='Ryu Hyunwoo',
      author_email='wowo0118@korea.ac.kr',
      description='Genome identification tool for Double-strand Break Repair',
      long_description=long_description,
      long_description_content_type='text/markdown',
      url='https://github.com/Chemical118/GDBr',
      packages=['gdbr'],
      package_dir={'gdbr' : 'gdbr/'},
      license='MIT',
      classifiers=['Development Status :: 4 - Beta',
                   'License :: OSI Approved :: MIT License'],
      install_requires=['RagTag', 'svim-asm', 'pyfaidx', 'vcfpy', 'biopython', 'requests'],
      python_requires='>=3.7',
      zip_safe=True,
      scripts=['gdbr/gdbr']
     )