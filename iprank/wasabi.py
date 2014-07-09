'''
Author: Alan medlar
'''
from commands import getstatusoutput as run
from sys import stderr, exit, argv

import webbrowser
import json


def wasabi(filename, name, url, uid):
    wid = _upload_to_wasabi(filename, name, url, uid)

    webbrowser.open_new_tab(_build_url(url, wid, uid))

def _build_url(url, wid, uid):
    return "%s/%s?share=%s" % (url, uid, wid)

def _upload_to_wasabi(filename, name, url, uid):
    curl_command = 'curl --silent \
                    -F "action=save" \
                    -F "writemode=new" \
                    -F "name=%s" \
                    -F "userid=%s" \
                    -F "file=@%s" %s' % (name, uid, filename, url)

    status,output = run(curl_command)

    if status != 0 :
        raise Exception("curl failed")

    tmp = json.loads(output)
    return tmp["id"]

if __name__ == "__main__":
	userid = argv[1]
	filename = argv[2]
	name = argv[3]
	url = "http://wasabi2.biocenter.helsinki.fi:8000"
	wasabi(filename, name, url, userid)
