
import sys, platform

def msg_exit(msg):
        print >> sys.stderr, msg+"\n"
        sys.exit(1)

if platform.system() == "Windows":
        msg_exit("Sorry, %s does not support Windows system so far."%package_name)
__version__ = "0.1.0"
package_name = "canopy"

bin_path = "bin"
global_path = '/usr/local/bin/' + package_name
local_path = '~/.local/bin/' + package_name
