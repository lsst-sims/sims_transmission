"""
MJDtools.py
--
Set of functions to convert between Modified Julian Date (MJD) and Gregorian date.
Also a function getSeason that returns the season (as an int, cf. function) for a
given MJD.
"""
import numpy as np
import datetime

MJDorigin = 678576

def fromMJD(mjd, numpy=None):
    """From MJD to Gregorian
    Input as float
    Output in the form [year, month, day, dayfraction]"""
    #mjd = np.asarray(mjd)
    #if len(mjd) == 1:
    if type(mjd) != np.ndarray:
        print mjd + MJDorigin
        date = datetime.date.fromordinal(int(mjd + MJDorigin))
        year, month, day = map(int, datetime.date.isoformat(date).split('-'))
    else:
        mjd_int = mjd.astype(int)
        date_arr = [datetime.date.fromordinal(mjd + MJDorigin) for mjd in mjd_int]
        ymd = [map(int, datetime.date.isoformat(date).split('-')) for date in date_arr]
        year, month, day = np.asarray(ymd, dtype='int').T
    tab = [year, month, day, np.mod(mjd,1.)]
    nptab = np.asarray(tab)
    return nptab 


def year2MJD(fyear):
    """From fractional year to MJD
    Input as float (or float array)
    Output as float (or float array)"""
    fyear_arr = np.asarray(fyear, dtype='float')
    year = np.asarray(fyear, dtype='int')
    days = (fyear_arr - year) * 365
    month = np.asarray(days // 30, dtype='int') + 1
    day = np.asarray(days % 30, dtype='int')
    dfrac = days - np.asarray(days, dtype='int')
    date_list = [datetime.date(y,m,d) for y,m,d in zip(year,month,day)]
    mjd_list = [d.toordinal() - MJDorigin for d in date_list]
    mjd_arr = np.asarray(mjd_list, dtype='int') + dfrac
    if mjd_arr.size==1:
        return float(mjd_arr)
    else:
        return mjd_arr

def toMJD(year, month, day, frac=None):
    """From date to MJD
    Input as year, month, day (or float arrays)
    Output as float (or float array)"""
    #print type(year), type(month), type(day)
    if type(year)==int or type(year)==np.float64:
        date_list = [datetime.date(int(year),int(month),int(day))]
    else:
        #print year
        year = np.asarray(year, dtype='int')
        month = np.asarray(month, dtype='int')
        day = np.asarray(day, dtype='int')
        #print year
        #print type(year)
        date_list = [datetime.date(y,m,d) for y,m,d in zip(year,month,day)]
    mjd_list = [d.toordinal() - MJDorigin for d in date_list]
    mjd_arr = np.asarray(mjd_list, dtype='int')
    if mjd_arr.size==1:
        if frac:
            return int(mjd_arr) + frac
        else:
            return int(mjd_arr)
    else:
        if frac:
            return mjd_arr + frac
        else:
            return mjd_arr

def getSeason(mjd, age=None, length=None):
    """Given MJD argument compute the season.
    Output :
    Season number (0,1,2,3) that corresponds to [North.hem. - South.hem.]
      0 : Winter  -  Summer
      1 : Fall    -  Spring
      2 : Summer  -  Winter
      3 : Spring  -  Fall
    Additional arguments:
    - Age : Days since the beginning of the season
    - Length : Season length
    """
    cur_year = fromMJD(mjd)[0]
    tab = np.ones_like(cur_year)
    if (mjd - toMJD(cur_year,tab*12,tab*21)) > 0 :
        cur_year += 1
    s0 = toMJD(cur_year-1,tab*12,tab*21)
    s1 = toMJD(cur_year,tab*3,tab*21)
    s2 = toMJD(cur_year,tab*6,tab*21)
    s3 = toMJD(cur_year,tab*9,tab*21)
    s4 = toMJD(cur_year,tab*12,tab*21)
    if (mjd-s4)>0:
        ssn = [0, mjd-s4-1]   #   Winter (Northern hem.) Summer (South)
    elif (mjd-s3)>0:
        ssn = [3, mjd-s3-1]   #   Spring (N) Fall (S)
    elif (mjd-s2)>0:
        ssn = [2, mjd-s2-1]   #   Summer (N) Winter (S)
    elif (mjd-s1)>0:
        ssn = [1, mjd-s1-1]   #   Fall (N)   Spring (S)
    else:
        ssn = [0, mjd-s0-1]   #   Winter (N) Summer (S)
    sld  = [s1-s0, s2-s1, s3-s2, s4-s3]
    ssn.append(sld[ssn[0]])
    if age:
        if length:
            return np.array(ssn)
        else:
            return np.array(ssn[:-1])
    elif length:
        return np.array([ssn[0], ssn[2]])
    else:
        return int(ssn[0])


def MJDtoindex(mjd, seas=None):
    """Get the day number in a normal or seasonal year
    given an mjd
    (seasonal year starts on the 21st of december)
    """
    started_year = fromMJD(mjd)[0]
    if seas:
        if (mjd - toMJD(started_year,12,21)) < 0 :
            started_year -= 1
        id_start = int(mjd - toMJD(started_year,12,21))
    else:
        id_start = int(mjd - toMJD(started_year,1,1))
    print "MJDtoindex:",mjd,  id_start    
    return id_start

