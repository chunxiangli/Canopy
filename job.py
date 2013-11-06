import os, sys, random, traceback, time
from StringIO import StringIO
from Queue import Queue
from threading import Thread, Event, Lock
from subprocess import Popen, PIPE
from logger import get_logger
from file_manage import open_file

_LOG = get_logger(__name__)

class JobQueue(Queue):
	def put(self, job):
		#_LOG.info("Job added:%s"%" ".join(job._command))
		Queue.put(self, job)

jobQueue = JobQueue()

class JobBase(object):
	def __init__(self, **kwargs):
		self.name = kwargs.get("description","whole")
		self._result = None	
		self.run_time = 0

		if "description" in kwargs:
			del kwargs["description"]

		self.id = kwargs.get("id", "")
		if "id" in kwargs:
			del kwargs["id"]

		self._kwargs = kwargs
		self._cwd = self._kwargs.get('cwd', os.curdir)

	@property
	def result(self):
		return self._result

	def start(self):
		pass

	def wait(self):
		pass
	
	def kill(self):
		pass

	def check_status(self):
		pass

class FakeJob(JobBase):
	def __init__(self, result, **kwargs):
		JobBase.__init__(self, **kwargs)
		self._result = result

	def start(self):
		self.run_time += 1

	def get_result(self):
		return self._result

class Job(JobBase):
	def __init__(self, command, post_processor, **kwargs):
		JobBase.__init__(self, **kwargs)
		self._command = command
		self._process = None
		self._return_code = None
		self._err = None
		self._stdout_stream = None
		self._stderr_stream = None
		self._post_processor = post_processor
		self._event_list = Event()
		self._finished_event = Event()
		self._waiting = False
		self._waiting_lock = Lock()

	@property
	def cwd(self):
		return self._cwd

	def start(self):
		self.run_time += 1

		assert self._process is None, "The process has existed."

		try:
			stdout_file = self._kwargs.get('stdout', None)
			stderr_file = self._kwargs.get('stderr', None)
			if stdout_file is not None:
				self._stdout_stream = open_file(stdout_file, 'w')
			else:
				self._stdout_stream = open_file(os.path.join(self._cwd, 'stdout.txt'), 'w')


			if stderr_file is not None:
				self._stderr_stream = open_file(stderr_file, 'w')
                        else:   
                                self._stderr_stream = open_file(os.path.join(self._cwd, 'stderr.txt'), 'w')

			self._kwargs['stdout'] = self._stdout_stream
			self._kwargs['stderr'] = self._stderr_stream 
			self._process = Popen(self._command, stdin = PIPE, **self._kwargs)
		except Exception as e:
			self._err = RuntimeError("Job Failed to start:%s. Erroinfo:%s"%(" ".join(self._command), str(e)))
			raise self._err
		finally:
			self._event_list.set()

	def wait(self):
		self._waiting_lock.acquire()	
		if self._waiting:
			self._waiting_lock.release()
			self.check_status()

			self._finished_event.wait()
			self.check_status()
		else:
			self._waiting = True
			self._waiting_lock.release()
			try:
				self.check_status()

                		if self._return_code is None:
                    			_LOG.debug('Waiting for process set up:%s'%self.name)
                    			self._event_list.wait()
                    			_LOG.debug('Job started:%s'%self.name)
					self.check_status()

                    		err_msg = []
                    		err_msg.append("Job failed: %s."%" ".join(self._command))

                    		try:
                        		self._return_code = self._process.wait()
					if self._return_code:
						err_msg.append("At:run")
                            			error = self.read_stderr()
                            			if error:
                                			err_msg.append(error)
                            			self._err = Exception("".join(err_msg))
                            			raise self._err
                    		except Exception as e:
					if self._err is None:
						err_msg.append("ErrorInfo:%s"%str(e))
						self._err = Exception("".join(err_msg))
                        		raise self._err
				finally:
					self._process.stdin.close()
                        		self._stdout_stream.close()
                        		self._stderr_stream.close()


                    		try:
                        		self._result = self._post_processor()
                    		except Exception, e:
                        		error = self.read_stderr()
					err_msg.append("At:post process.")

                        		if not error:
                            			error = str(e)

                        		err_msg.append(error)
                        		self._err = Exception("".join(err_msg))
				self.check_status()
            		finally:
                		self._finished_event.set()

            	return self._return_code

	def read_stderr(self):
		if os.path.exists(self._stderr_stream.name):
			stderr_stream = open(self._stderr_stream.name, 'r')
			error = stderr_stream.read(-1)	
			stderr_stream.close()
			return error
		else:
			return None

	def check_status(self):
		if self._err is not None:
            		raise self._err

	def get_result(self):
		self.check_status()
        	if self._result is None:
            		self.wait()

        	return self._result

    	def kill(self):
        	#if self._result is None:
		if self._process is not None:
            		self._process.kill()
                                  

class Worker(Thread):
	id = 0
	def __init__(self):
		Thread.__init__(self)
		self._exitFlag = False
		Worker.id += 1
		self.id = Worker.id

	def run(self):
		while True:
			if not jobQueue.empty():
				job = jobQueue.get()
				try:
					_LOG.info("Worker %d: Job start: %s."%(self.id, job.name))
					job.start()
				except Exception as e:
					#TODO:whether need to rerun job
					if job.run_time < 1:
						_LOG.info("reput job into queue")
						jobQueue.put(job)
					else:
						self.log_error(str(e))
						raise e
				else:
					try:
						job.get_result()
					except Exception as e:
						if job.run_time < 1:
							jobQueue.put(job)
						else:
							self.log_error(str(e))
							raise e
					else:
						_LOG.info("Worker %d: Job finished: %s."%(self.id, job.name))
						jobQueue.task_done()
	def exit(self):
		self._exitFlag = True

	def log_error(self, msg):
		msg = "Worker %d: %s."%(self.id, msg)
                _LOG.error(msg)

class MainWorker(object):
	def __init__(self, max_num=1):
		self._workers_lock = Lock()
		self._worker_list = []
		self._max_num_workers = max_num

	def set_max_num_workers(self, num):
		self._workers_lock.acquire()
		cur_num = len(self._worker_list)

		if num < self._max_num_workers and num < cur_num:
			for i in range(num):
				self._worker_list[i].exit()	

		self._max_num_workers = num
		self._workers_lock.release()

	def start_workers(self, num=1):
		self._workers_lock.acquire()
		cur_num = len(self._worker_list)
		new_num = cur_num + num

		if new_num > self._max_num_workers:
			new_num = self._max_num_workers
			_LOG.warn("%d workers are running. Cannot start another %d workers. Maximum number %d workers will be running.."%(cur_num, num, self._max_num_workers))

		for i in range(cur_num, new_num):
			w = Worker()
			self._worker_list.append(w)
			'''
			Program will terminated when no active non-deamon thred. Don't call sys.exit 
			when you still need operate on modules and global variables in the global namespace.
			'''
			w.setDaemon(True)
			w.start()
		self._workers_lock.release()

mainWorker = MainWorker()
