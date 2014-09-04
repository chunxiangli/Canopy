from distutils.core import setup

setup(name='canopy',
      version='1.0',
      description='Co-estimation alignment and phylogeny in Python',
      author='Chunxiang Li',
      author_email='chunxiang.li@helsinki.fi',
      platforms=['*nix'],
      requires=['Biopython(>=1.58)', 'Dendropy(>=3.10.0)'],
      packages=['canopy'],
      scripts=['scripts/canopy']
     )
