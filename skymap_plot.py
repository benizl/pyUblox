#!/usr/bin/env python

import sys
import numpy, math
import scipy.sparse
import pylab as plt

import ephemeris, satPosition, util

from optparse import OptionParser

import util, ublox

parser = OptionParser("skymap_plot.py [options]")
parser.add_option("--ubx-log", help="uBlox log file from which Ephemerides are loaded") 

(opts, args) = parser.parse_args()

t_first = 0
t_last = 0
t_wrap = 0

sat_el = scipy.sparse.lil_matrix((200*60*60, 140))
sat_az = scipy.sparse.lil_matrix((200*60*60, 140))

print("Parsing UBX")
dev = ublox.UBlox(opts.ubx_log)

t = 0
eph = {}

azs = [[] for i in range(33)]
els = [[] for i in range(33)]

ourpos = None

while True:
    '''process the ublox messages, extracting the ones we need for the sat position'''
    msg = dev.receive_message()
    if msg is None:
        break
    if msg.name() in [ 'AID_EPH', 'RXM_EPH' ]:
        msg.unpack()
        eph[msg.svid] = ephemeris.EphemerisData(msg)
    elif msg.name() == 'RXM_RAW':
        if ourpos is None:
            continue

        msg.unpack()
        t = msg.iTOW * 0.001

        for s in msg.recs:
            svid = s.sv
            if svid not in eph or svid > 32:
                continue

            satpos = satPosition.satPosition_raw(eph[svid], svid, t)
            if satpos is None:
                continue

            az, el = satPosition.calculateAzimuthElevation_raw(satpos, ourpos)

            if el > 5:
                azs[svid].append(math.radians(az))
                els[svid].append(math.cos(math.radians(el)))
    elif msg.name() in ['NAV_SOL', 'NAV_POSECEF']:
        msg.unpack()
        ourpos = util.PosVector(msg.ecefX * 0.01, msg.ecefY * 0.01, msg.ecefZ * 0.01)


plt.figure()
plt.suptitle("Skyviews")
ax = plt.subplot(1,1,1,polar=True)
ax.set_rmax(1.0)
ax.set_theta_offset(math.pi / 2)
ax.set_theta_direction(-1)

for sat in range(32):
    if len(azs[sat]) > 0:
        ax.plot(azs[sat], els[sat])

plt.show()


