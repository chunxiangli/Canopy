from distutils.core import setup

setup(name='iprank',
      version='1.0',
      description='Iterative alignment with PRANK',
      author='Chunxiang Li',
      author_email='chunxiang.li@helsinki.fi',
      platforms=['*nix'],
      requires=['Biopython(>=1.58)', 'Dendropy(>=3.10.0)'],
      packages=['iprank'],
      scripts=['scripts/iPRANK']
     )
