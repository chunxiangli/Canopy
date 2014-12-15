import sys, platform

def msg_exit(msg):
        print >> sys.stderr, msg+"\n"
        sys.exit(1)

__version__ = "0.1.1"
package_name = "canopy"

if platform.system() == "Windows" or platform.system() == "Darwin":
        msg_exit("Sorry, %s only support Unix-like system so far."%package_name)

bin_path = "bin"
global_path = '/usr/local/bin/' + package_name
local_path = '~/.local/bin/' + package_name

