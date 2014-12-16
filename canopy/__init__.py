import sys, platform
import site

def msg_exit(msg):
        print >> sys.stderr, msg+"\n"
        sys.exit(1)

__version__ = "0.1.02"
package_name = "canopy"

if platform.system() == "Windows":
        msg_exit("Sorry, %s does not support Windows system so far."%package_name)

bin_path = "bin"
if platform.system() == "Darwin":
                bin_path = "mac_bin"

global_path = '/usr/local/bin/' + package_name
local_path = '%s/bin/%s'%(site.USER_BASE, package_name)
