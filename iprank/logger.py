import os, sys, logging

class Messenger(object):
    """
    Wraps reporting of messages, progress, warnings and errors to user.
    Singleton (instantiated below).
    """

    def __init__(self):
        self.err_log_streams = [sys.stderr]
        self.out_log_streams = [sys.stdout]

    def _format_message(self, msg, msg_type):
        return "%s: %s\n" % (msg_type, msg)

    def _write_to_streams(self, streams, msg, flush=True):
        for s in streams:
            s.write(msg)
            if flush:
                s.flush()

    def send_error(self, msg):
        msg = self._format_message(msg, "ERROR")
	self._write_to_streams(self.err_log_streams, msg)

    def send_warning(self, msg):
        msg = self._format_message(msg, "WARNING")
	self._write_to_streams(self.err_log_streams, msg)

    def send_info(self, msg):
        msg = self._format_message(msg, "INFO")
	self._write_to_streams(self.out_log_streams, msg)

MESSENGER = Messenger()


_LOGGING_LEVEL_EVAR = "LOGGING_LEVEL"
_LOGGING_FORMAT_EVAR = "LOGGING_FORMAT"


def get_logging_level():
    """
	Check environment for LOGGING_LEVEL and return a logging level
    integer.
    """
    if _LOGGING_LEVEL_EVAR in os.environ:
        if os.environ[_LOGGING_LEVEL_EVAR].upper() == "NOTSET":
            level = logging.NOTSET
        elif os.environ[_LOGGING_LEVEL_EVAR].upper() == "DEBUG":
            level = logging.DEBUG
        elif os.environ[_LOGGING_LEVEL_EVAR].upper() == "INFO":
            level = logging.INFO
        elif os.environ[_LOGGING_LEVEL_EVAR].upper() == "WARNING":
            level = logging.WARNING
        elif os.environ[_LOGGING_LEVEL_EVAR].upper() == "ERROR":
            level = logging.ERROR
        elif os.environ[_LOGGING_LEVEL_EVAR].upper() == "CRITICAL":
            level = logging.CRITICAL
        else:
            level = logging.NOTSET
    else:
        level = logging.NOTSET

    return level

def get_logger(name="iterative_prank"):
    """
    Returns a logger with name set as given, and configured
    to the level given by the environment variable _LOGGING_LEVEL_ENVAR.
    """
    logger_set = False
    logger = logging.getLogger(name)

    if not logger_set:
        level = get_logging_level()

        rich_formatter = logging.Formatter("[%(asctime)s] %(filename)s (line %(lineno)d): [%(levelname) 8s]: %(message)s")
        simple_formatter = logging.Formatter("[%(levelname) 8s]: %(message)s")
        default_formatter = logging.Formatter("%(name) 11s:[%(levelname)s]:%(message)s")
        logging_formatter = default_formatter

        if _LOGGING_FORMAT_EVAR in os.environ:
            if os.environ[_LOGGING_FORMAT_EVAR].upper() == "RICH":
                logging_formatter = rich_formatter
            elif os.environ[_LOGGING_FORMAT_EVAR].upper() == "SIMPLE":
                logging_formatter = simple_formatter
            elif os.environ[_LOGGING_FORMAT_EVAR].upper() == "NONE":
                logging_formatter = None
            else:
                logging_formatter = default_formatter
        else:
            logging_formatter = default_formatter

        if logging_formatter is not None:
            logging_formatter.datefmt = '%H:%M:%S'

        logger.setLevel(level)
        ch = logging.StreamHandler()
        ch.setLevel(level)
        ch.setFormatter(logging_formatter)
        logger.addHandler(ch)

    return logger
