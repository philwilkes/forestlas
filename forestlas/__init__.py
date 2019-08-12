# http://www.diveintopython3.net/case-study-porting-chardet-to-python-3.html#multifile-modules

def detect(aBuf):
    from . import universaldetector
    u = universaldetector.UniversalDetector()
    u.reset()
    u.feed(aBuf)
    u.close()
    return u.result