#!/usr/bin/env python
import requests
import sys
import os.path
import subprocess

date = sys.argv[1]
fields = set(int(arg) for arg in sys.argv[2:])
filt = 'R'

r = requests.post(
    'http://ptfdepot.ipac.caltech.edu/cgi-bin/getNightlyRealTimeImages.cgi',
    data={'nightdate': date, 'outputformat': 'JSON', 'action': 'Query Database'}).json()

for row in r:
    if row['ptffield'] in fields and row['filter'] == 'R':
        imgname = os.path.basename(row['filename'])
        catname = imgname.replace('_flattened.fits', '_sex.cat')
        for filename in [imgname, catname]:
            url = 'http://ptfdepot.ipac.caltech.edu/ptfdata/RealTimeProds/' + date.replace('-', '') + '/' + filename
            subprocess.call(['wget', '-N', url])
