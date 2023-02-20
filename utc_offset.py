# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 22:00:16 2022

@author: Zoletnik
"""

from datetime import datetime, timedelta, timezone
import pytz

def UTC_offset(UTC_offset_minutes=None,date=None,year=None,month=None,day=None):
    """
    Returns the UTC offset.

    Parameters
    ----------
    UTC_offset_minutes : int, optional
        The UTC offset in minutes. If none of the other paramteres is set this will be returned. The default is None.
    date : string, optional
        The date YYYYMMDD. The default is None. If this is set it is used instead of year,month,day
    year : string or int, optional
        The year. The default is None.
    month : string or int, optional
        The month. The default is None.
    day : string or int, optional
        The day. The default is None.

    Raises
    ------
    ValueError
        Bad date format.

    Returns
    -------
    int
        The UTC offset in minutes.

    """
    if ((date is None) and (year is None) and (month is None) and (day is None)):
        return UTC_offset_minutes
    if (type(date) is str):
        if (len(date) != 8):
            raise ValueError("date argument in UTC_offset should have format YYYYMMDD.")
        _year = int(date[:4])
        _month = int(date[4:6])
        _day = int(date[6:8])
    else:
        if ((year is None) or (month is None) or (day is None)):
            raise ValueError('If date is not set string year,month, day should be set.')
        _year = int(year)
        _month = int(month)
        _day = int(day)
    timezone = pytz.timezone("CET")
    return timezone.localize(datetime(_year,_month,_day,12,0,0)).utcoffset() // timedelta(minutes=1)