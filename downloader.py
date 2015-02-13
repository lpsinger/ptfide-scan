#!/usr/bin/env python
import requests
import sys
import os.path
import shutil

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
            print filename
            url = 'http://ptfdepot.ipac.caltech.edu/ptfdata/RealTimeProds/' + date.replace('-', '') + '/' + filename
            r = requests.get(url, stream=True)
            r.raise_for_status()
            r.raw.decode_content = True
            with open(filename, 'wb') as f:
                shutil.copyfileobj(r.raw, f)
