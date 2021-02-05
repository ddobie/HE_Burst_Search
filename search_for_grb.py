""" 
Given a specific position and time and search radius
and search window (in time),
check IPN, Fermi, and Swift to see whether any GRBs
were detected.

Note that this assumes that the year is 20XX.
"""

import numpy as np
import requests
import lxml.html as lh
import pandas as pd
import re
import subprocess
import os
from datetime import datetime
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u

def get_searchstr(t):
    """ useful function for the IPN catalog search """
    yy = str(t.datetime.year)[2:]
    mm = t.strftime('%b').upper()
    dd = str(t.datetime.day).zfill(2)
    searchstr = ' '.join([dd,mm,yy])
    return searchstr

######################################

def ipn(ra,dec,start,end):
    """ Check the online IPN catalog
    Note that this catalog is not updated in real time

    This code will only work for bursts in 2000 onward.
    Would need to modify for years with 19XX.
    """
    c = SkyCoord(ra,dec,unit='deg')

    print("CONDUCTING SEARCH OF IPN CATALOG")

    # convert to JD
    window = [Time(start, format='datetime'), Time(end, format='datetime')]

    # Make sure that each day is represented
    window_grid = Time(
            np.arange(window[0].jd, window[-1].jd+1, 0.5), format='jd')
    searchstr = [get_searchstr(t) for t in window_grid]
    searchstr = np.unique(np.array(searchstr))

    # Pull out the relevant lines
    fpath = "http://www.ssl.berkeley.edu/ipn3/masterli.txt"
    fname = fpath.split('/')[-1]
    subprocess.call(['wget', '-O', 'masterli.txt', 'http://www.ssl.berkeley.edu/ipn3/masterli.txt'])
    lines = np.array(open(fname, "r").readlines())
    header = lines[np.array([' DOY TIME ' in l for l in lines])][0]
    keep = np.array([l[7:16] in searchstr for l in lines])

    # Now, check each one to see if the time is correct
    final_set = []
    for l in lines[keep]:
        filtered_l = [i for i in l.split(" ") if i]
        burst_dd = filtered_l[0].split('.')[1]
        burst_mm = str(strptime(filtered_l[1], '%b').tm_mon).zfill(2)
        burst_yy = str(filtered_l[2]).zfill(2)
        burst_time = filtered_l[4]
        burst_datetime = Time(
                '20%s-%s-%sT%s' %(burst_yy,burst_mm,burst_dd,burst_time), 
                format='isot')
        if np.logical_and(
                burst_datetime >= window[0], burst_datetime <= window[-1]):
            final_set.append(l)

    print("There are %s bursts in the %s-day window" %(
        len(final_set), window[-1]-window[0]))

    # Check which spacecraft observed these bursts
    for l in final_set:
        det_by = np.array(
                [header[i.start():i.end()] for i in re.finditer('YES', l)])
        print(det_by)

######################################

def fermi(ra,dec,start,end):
    """ Check the online Fermi trigger catalog 

    Under the 'trigger_type' entry, GRBs are under 'GRB' and 
    the 'reliability' entry gives the classification probability. 
    Fermi usually does >80% GRB for their automatic processing threshold
    """
    c = SkyCoord(ra,dec,unit='deg')
    window = [Time(start, format='datetime'), Time(end, format='datetime')]

    print("\n")
    print("CONDUCTING SEARCH OF FERMI CATALOG")

    www = 'https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3query.pl'
    data = {}
    data['tablehead'] = 'name=heasarc_fermigtrig&description=Fermi GBM Trigger Catalog&url=http://heasarc.gsfc.nasa.gov/W3Browse/fermi/fermigtrig.html&archive=Y&radius=180&mission=FERMI&priority=1&tabletype=Object'
    data['Time'] = '%s .. %s' %(window[0].mjd, window[-1].mjd)
    data['displaymode'] = 'PureTextDisplay'
    data['varon'] = 'name,ra,dec,trigger_time,error_radius'
    r = requests.post(url = www, data = data)
    out = np.array([i for i in r.text.split('\n') if i])
    header = [i for i in out[2].split('|')]
    ncands = len(out[3:])
    if ncands == 0:
        print("No GRBs in Fermi")
    else:
        print("Found %s in Fermi" %ncands)
        for ii in np.arange(ncands):
            cand = out[3+ii]
            vals = [i for i in cand.split('|')]
            grbname = vals[1]
            grbra = vals[2]
            grbdec = vals[3]
            grbtime = vals[4]
            grbepos = vals[5].strip()
            print("%s with RA=%s, Dec=%s on t=%s with %s deg uncertainty"%(
                grbname,grbra,grbdec,grbtime,grbepos))
            c2 = SkyCoord('%s %s' %(grbra,grbdec), unit=(u.hourangle, u.deg))
            dist = c2.separation(c).degree
            print("The burst is %s deg away from the source" %dist)


def fermi_subthreshold(ra,dec,start,end):
    """ 
    check the fermi subthreshold notices
    https://gcn.gsfc.nasa.gov/fermi_gbm_subthresh_archive.html

    the list only starts in 2017. to search it prior to that date,
    you would have to do something different.

    the html scraper borrows heavily from
    https://towardsdatascience.com/web-scraping-html-tables-with-python-c9baba21059
    """

    window = [Time(start, format='datetime'), Time(end, format='datetime')]
    c = SkyCoord(ra,dec,unit='deg')

    print('\n')
    print("CONDUCTING SEARCH OF FERMI SUBTHRESHOLD CATALOG")
    url = "https://gcn.gsfc.nasa.gov/fermi_gbm_subthresh_archive.html"
    page = requests.get(url)
    doc = lh.fromstring(page.content)
    tr_elements = doc.xpath('//tr')

    # the first element is notes
    # the second element is the header
    col = []
    for ii,t in enumerate(tr_elements[1]):
        name = t.text_content()
        col.append((name,[]))

    # the third element onwards are rows of the table
    for j in range(2, len(tr_elements)):
        T = tr_elements[j]
        i = 0
        for t in T.iterchildren():
            data = t.text_content()
            col[i][1].append(data)
            i+=1

    # create dataframe
    Dict={title:column for (title,column) in col}
    df=pd.DataFrame(Dict)

    # Now, go through the table and check whether each burst
    # is consistent with the time interval.
    # If so, calculate the RA and Dec offset
    for index,row in df.iterrows():
        date = '20' + str('-'.join(row['Date'].split('/')))
        time = row['Time UT']
        datetime = Time('%sT%s' %(date,time), format='isot')
        if np.logical_and(datetime > window[0], datetime < window[-1]):
            if row['Rel'] != '2':
                ra = row['RA(J2000)[deg]']
                dec = row['Dec(J2000)[deg]']
                epos = row['Error[deg]']
                print("Subthreshold burst at %s pos %s, %s with error %s" %
                        (datetime, ra, dec, epos))
                c2  = SkyCoord(ra, dec, unit='deg')
                dist = c.separation(c2).degree
                print("Distance is %s deg from source" %dist)
                return ""
    print("No Fermi subthreshold bursts found")

######################################

def swift(ra,dec,start,end):
    """ 
    check the Swift GRB table
    https://swift.gsfc.nasa.gov/archive/grb_table/

    sometimes, the burst time is listed as nan. that's when a burst was found in subsequent
    ground analysis. for now I am ignoring this case, since typically bursts were found
    by other spacecraft.
    """
    c = SkyCoord(ra,dec,unit='deg')

    print("\n")
    print("CONDUCTING SEARCH OF SWIFT/BAT CATALOG")

    www = 'https://swift.gsfc.nasa.gov/archive/grb_table/table.php?'
    info = "obs=Swift&year=All+Years&grb_time=1&bat_ra=1&bat_dec=1&bat_err_radius=1"
    df = pd.read_html(www+info)[0]
    
    names = np.array([val[0][0:6] for val in df['GRB'].values])
    ttime = np.array([val[0] for val in df['Time[UT]'].values])
    grbra = np.array([val[0] for val in df['BAT RA(J2000)'].values])
    grbdec = np.array([val[0] for val in df['BAT Dec(J2000)'].values])
    grbepos = np.array([val[0] for val in df['BAT 90%Error Radius[arcmin]'].values])

    # ignore the bursts where the time is n/a or the pos is n/a
    ignore = np.logical_or(ttime=='nan', grbra=='nan')
    names = names[~ignore]
    ttime = [val[0:8] for val in ttime[~ignore]] # get rid of decimal places
    grbt_raw = np.array([''.join(map(str, i)) for i in zip(names,ttime)])
    grbt = np.array([datetime.strptime(val, '%y%m%d%H:%M:%S') for val in grbt_raw])
    grbra = [val[0:6] for val in grbra[~ignore]]
    grbdec = [val[0:6] for val in grbdec[~ignore]]
    grbepos = grbepos[~ignore]

    keep = np.logical_and(grbt >= start, grbt <= end)
    if sum(keep) == 0:
        print("No Swift/BAT bursts during this window")
    else:
        print("%s Swift/BAT bursts during this window" %sum(keep))
        for i in np.where(keep)[0]:
            print("Position is %s, %s with error %s deg" %(
                grbra[i],grbdec[i],grbepos[i]))
            c2 = SkyCoord(grbra[i], grbdec[i], unit='deg')
            print("Separation is %s deg" %c.separation(c2).degree)


def run(ra,dec,start,end):
    c = SkyCoord(ra, dec, unit='deg')
    ipn(ra,dec,start,end)
    fermi(ra,dec,start,end)
    fermi_subthreshold(ra,dec,start,end)
    swift(ra,dec,start,end)


if __name__=="__main__":
    start = datetime.fromisoformat('2020-10-08T00:00:00')
    end = datetime.fromisoformat('2020-10-13T00:00:00')
    ra = 335.008456
    dec = -2.840352
    run(ra,dec,start,end)
