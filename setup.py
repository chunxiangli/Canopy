import glob, sys, os
from distutils.core import setup
import canopy

target_dir = canopy.global_path
make_symlink = False
user_path = os.path.expanduser('~/bin')
symlink_path = os.path.join(user_path, canopy.package_name)
if '--user' in sys.argv:
        if not os.path.isdir(user_path):
                canopy.msg_exit('%s does not exist! Please make sure it included in the environment variable "PATH".'%user_path) 
        if os.path.exists(symlink_path):
                canopy.msg_exit('%s existed!\nYou should remove it and then try again.'%symlink_path)     

	make_symlink = True
	target_dir = os.path.expanduser(canopy.local_path)

target_dir += canopy.bin_path
extend_files = {target_dir:[]}
source_dir = os.path.join(canopy.package_name, canopy.bin_path)

def walk_dir(sdir, tdir):
        for f in glob.glob(os.path.join(sdir, '*')):
                if os.path.isfile(f):
                        extend_files.setdefault(tdir, []).append(f)
                else:
                        subdir = os.path.relpath(f, sdir)
                        ntdir = os.path.join(tdir, subdir)
                        walk_dir(f, ntdir)
walk_dir(source_dir, target_dir)

setup(name=canopy.package_name,
      version=canopy.__version__,
      description='Co-estimation alignment and phylogeny in Python',
      author='Chunxiang Li',
      url='http://wasabiapp.org/software/canopy',
      license='GNU GENERAL PUBLIC LICENSE Version 3(http://www.gnu.org/licenses/)',
      author_email='chunxiang.li@helsinki.fi',
      platforms=['*nix', 'Darwin'],
      requires=['Biopython(>=1.58)', 'Dendropy(>=3.10.0)'],
      packages=[canopy.package_name],
      data_files= sorted(extend_files.items()),
      scripts=['scripts/'+canopy.package_name]
     )

if make_symlink:
        os.symlink(os.path.expanduser(canopy.local_path), symlink_path)	
