import os, os.path
from contextlib import contextmanager
import glob

@contextmanager
def open(path):
	file_path = os.path.join(path, 'assembly.log')
	if os.path.exists(os.path.join(path, 'assembly.log')):
		if os.stat(file_path).st_size > 102400: #100kb
			os.rename(file_path, os.path.join(path, 'assembly.{}.log'.format(
				len(glob.glob(os.path.join(path, 'assembly*.log'))))))
	file_handle = __builtins__['open'](os.path.join(path, 'assembly.log'), 'a')
	yield file_handle
	file_handle.flush()
	file_handle.close()