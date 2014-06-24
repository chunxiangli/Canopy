import sys, os, shutil, tarfile
from threading import Lock
from iprank.logger import get_logger

_LOG = get_logger(__name__)


def remove_files(file_list):
	_LOG.debug("remove files:%s"%file_list)
	try:
		for f in file_list:
			os.remove(f)
	except(IOError):
		raise Exception("remove file %s failed."%f)

def remove_files_from_dir(work_dir, keyword):
	assert os.path.isdir(work_dir), "The %s isn't directory."%work_dir
	assert os.path.exists(work_dir), "The directory %s doesn't exist."%work_dir

	delete_files = [ os.path.join(work_dir, f) for f in os.listdir(work_dir) if keyword in f ]
	if not delete_files:
		raise IOError("No file contains %s"%keyword)
	else:
		remove_files(delete_files)
		
def open_file(file_path, mode=None):
	full_file_path = os.path.abspath(file_path)
	d = os.path.dirname(full_file_path)
	if not os.path.exists(d):
		os.makedirs(d)
	return open(full_file_path, mode)	
		
def copy_files(work_dir, keyword, new_keyword, target_dir=None):
	assert os.path.exists(work_dir), "The directory %s doesn't exist."%work_dir

	if target_dir is None:
		target_dir = work_dir
        source_files = filter(lambda x: keyword in x, os.listdir(work_dir))

        assert len(source_files), "No file name contains  %s."%keyword

        for f in source_files:
        	prefix, suffix = f.split(keyword)
         	new_fname = os.path.join(target_dir, new_keyword+suffix)
                shutil.copy2(os.path.join(work_dir, f), new_fname)

def makeArchive(target, filter_arr=None):
	full_path = os.path.realpath(target)
	if not os.path.exists(full_path):
		raise OSError("%s directry not exists."%full_path)

        base_name = os.path.basename(full_path)
        base_path = full_path[:-len(base_name)-1]
	origin_dir = os.getcwd()
	os.chdir(base_path)
        len_base = len(base_path)
	target_name = os.path.join(base_path, base_name+".tgz")
        target_file = tarfile.open(target_name, 'w:gz')
	if os.path.isfile(full_path):
		target_file.add(base_name)	
	else:
        	for dirpath, dirs, fs in os.walk(full_path):
         		relative_path = dirpath[len_base+1:]
                	if len(fs) > 0:
                		for f in fs:
                        		if filter_arr is None or not (os.path.join(dirpath, f) in filter_arr):
                                		target_file.add(os.path.join(relative_path, f))

        target_file.close()
	os.chdir(origin_dir)
	return target_name


class TempFileManager(object):
	'''
		All created temp files should be deleted by client code calling the same TempFileManager instance.
	'''
	def __init__(self, top_level_temp_real=None):
		self._directories_created = set()
		self._directories_created_lock = Lock()
		self._top_level_temp_real = top_level_temp_real
	
	def get_directories_created(self):
		return self._directories_created

	directories_created = property(get_directories_created)

	def _is_already_created(self, real_path):
		self._directories_created_lock.acquire()
		b = real_path in self._directories_created
		slef._directories_created_loc.release()
		return b

	def create_top_level_temp(self, parent, temp_name="temp"):
		assert self._top_level_temp_real is None, "Top level temporary directory '%s' exists."%self._top_level_temp_real

		parent_real = os.path.realpath(parent)
		if not os.path.exists(parent_real):
			raise OSError("Parent path %s doesn't exists."%parent_real)

		if not os.path.isdir(parent_real):
			raise OSError("Parent path %s is not a directory."%parent_real)

		top_level_temp = os.path.join(parent_real, temp_name)
		os.makedirs(top_level_temp)
		self._top_level_temp_real = os.path.realpath(top_level_temp)
		self._directories_created.add(self._top_level_temp_real)

		return self._top_level_temp_real

	def _is_subdir(self, subdir_real):
		if self._top_level_temp_real is None:
			raise OSError("Top level temproray directory has not been created or set.")

		common_str = os.path.commonprefix([self._top_level_temp_real, subdir_real])
		return common_str == self._top_level_temp_real

	def create_subdir(self, dir_name):
		'''
			create a directory under the top level temporary directory.
		'''
		dir_real_path = os.path.join(self._top_level_temp_real, dir_name)
		if os.path.exists(dir_real_path):
			raise OSError("Subpath '%s' already exists."%dir_real_path)

		self._directories_created_lock.acquire()
		try:
			if dir_real_path in self._directories_created:
				raise OSError("Subpath '%s' is still in _directories_created set"%dir_real_path)

			os.makedirs(dir_real_path)
			self._directories_created.add(dir_real_path)
			return dir_real_path
		finally:
			self._directories_created_lock.release()

	def create_temp_subdir(self, parent, dir_name):
		parent_real = os.path.realpath(parent)
		if not os.path.exists(parent_real):
			raise OSError("Parent '%s' for temp_subdir doesn't existed."%parent_real)

		if not os.path.isdir(parent_real):
			raise OSError("Parent '%s' for temp_subdir is not a directory."%parent_real)

		if not self._is_subdir(parent_real):
			raise OSError("Parent '%s' for temp_subdir is not under top level temporary directory '%s'."%(parent_real, self._top_level_temp_real))

		self._directories_created_lock.acquire()
		try:
			subdir_path = os.path.join(parent_real, dir_name)
			os.makedirs(subdir_path)
			subdir_path_real = os.path.realpath(subdir_path)
			self._directories_created.add(subdir_path_real)
			return subdir_path_real

		finally:
			self._directories_created_lock.release()

	def remove_dirs(self, del_dirs):
		for del_dir in del_dirs:
			dir_real_path = os.path.realpath(del_dir)
			if os.path.exists(dir_real_path):
				shutil.rmtree(dir_real_path)

			self._directories_created_lock.acquire()
			try:
				if dir_real_path in self._directories_created:
					self._directories_created.remove(dir_real_path)
					if dir_real_path == self._top_level_temp_real:
						self._top_level_temp_real = None
				else:
					raise OSError("'%s' is not a registered temp directory created by this processing!"%dir_real_path)
					
			finally:
				self._directories_created_lock.release()

	def clear(self):
                if self._top_level_temp_real is not None:
                        self.remove_dirs([self._top_level_temp_real])

        def makeArchive(self, filter_arr=None):
		target = makeArchive(self._top_level_temp_real, filter_arr=filter_arr)
                self.clear()
		return target
