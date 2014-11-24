import glob, sys
import os
import os.path as path
from distutils.core import setup
import canopy

target_dir = canopy.global_path
make_symlink = False
if '--user' in sys.argv:
	make_symlink = True
	target_dir = path.expanduser(canopy.local_path)
target_dir += canopy.bin_path
extend_files = {target_dir:[]}
source_dir = path.join(canopy.package_name, canopy.bin_path)

for f in glob.glob(path.join(source_dir,'*')):
	if path.isfile(f):
		extend_files[target_dir].append(f)
	else:
		subdir = path.relpath(f, source_dir)
		for sf in glob.glob(path.join(f, '*')):
			if path.isfile(sf):
				extend_files.setdefault(path.join(target_dir, subdir), []).append(sf)
setup(name=canopy.package_name,
      version=canopy.__version__,
      description='Co-estimation alignment and phylogeny in Python',
      author='Chunxiang Li',
      author_email='chunxiang.li@helsinki.fi',
      platforms=['*nix'],
      requires=['Biopython(>=1.58)', 'Dendropy(>=3.10.0)'],
      packages=[canopy.package_name],
      data_files= sorted(extend_files.items()),
      scripts=['scripts/'+canopy.package_name]
     )

if make_symlink:
	os.symlink(path.expanduser(canopy.local_path), path.join(path.expanduser('~'), 'bin', canopy.package_name))	
