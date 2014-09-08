#import subprocess
import tempfile
from contextlib import contextmanager
import uuid
import os
import shutil

"""
	
	Tempfile overrides to use ramdisks.

"""

def mkdtemp(hint=2147483648): #2 gigabytes
	if (os.statvfs('/run/shm').f_bfree*os.statvfs('/run/shm').f_bsize) > hint:
		tempdir_name = os.path.join('/dev/shm', str(uuid.uuid4()))
		os.mkdir(tempdir_name)
		return tempdir_name
	return tempfile.mkdtemp()
	
@contextmanager
def TemporaryFile(mode='w+b'):
	temp_dir = mkdtemp()
	filehandle = open(os.path.join(temp_dir, str(uuid.uuid4())), mode)
	yield filehandle
	filehandle.close()
	shutil.rmtree(temp_dir)
	
NamedTemporaryFile = TemporaryFile