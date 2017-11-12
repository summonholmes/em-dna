#!/usr/bin/env python
# This file is only to make the file compatible on this system.  It can be ignored.
with open('main.py', 'rb+') as f:
    content = f.read()
    f.seek(0)
    f.write(content.replace(b'\r', b''))
    f.truncate()
