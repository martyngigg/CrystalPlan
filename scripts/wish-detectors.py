from __future__ import (absolute_import, division, print_function, unicode_literals)
from mantid.simpleapi import *
import numpy as np
from mantid import mtd

ws = mtd["37828_hexaferrite_vana"]
instrument = ws.getInstrument()

with open("/Users/samueljackson/Desktop/WISH_full.csv", 'w') as fhandle:
    fhandle.write("#Cyl\n")

    header = "#{:24}  {:10}  {:10}  {:20}  {:20}  {:20}  {:20}  {:20}\n".format("name", "x", "y", "z", "radius(cm)", "height(cm)", "start_angle", "end_angle")
    fhandle.write(header)

    for i in range(1, 11):
        name = str("WISH/panel{:02}/WISHpanel{:02}".format(i, i))
        panel_name = str("Panel {}".format(i))
        if i == 6:
            name = str("WISH/panel{:02}/WISHpanel01".format(i, i))

        sample = instrument.getComponentByName(str('nickel-holder'))
        panel = instrument.getComponentByName(name)

        pixel1 = instrument.getComponentByName(str(name + "/tube001/pixel0001"))
        pixel2 = instrument.getComponentByName(str(name + "/tube001/pixel0512"))


        height = (pixel2.getPos() - pixel1.getPos()).norm()

        pixel3 = instrument.getComponentByName(str(name + "/tube152/pixel0001"))

        p1 = np.array(pixel1.getPos())
        p1 /= np.linalg.norm(p1)

        p2 = np.array(pixel3.getPos())
        p2 /= np.linalg.norm(p2)

        def angle_between(b):
            angle = np.arctan2(b[0], b[2])
            return angle


        radius = panel.getDistance(sample)

        start_angle = angle_between(p1)
        end_angle = angle_between(p2)

        if start_angle > end_angle:
            start_angle, end_angle = end_angle, start_angle

        #start_angle += np.pi/2.
        #end_angle += np.pi/2.
        #sstart_angle, end_angle = np.degrees(start_angle), np.radians(end_angle)
        print(str("{} start: {} end: {}".format(panel_name, np.degrees(start_angle), np.degrees(end_angle))))
        radius *= 100
        height *= 100

        line = "{:25}, {:10}, {:10}, {:20}, {:20}, {:20}, {:20}, {:20}".format(panel_name, 0, -height/2., 0, radius, height, start_angle, end_angle)
        fhandle.write(line)
        fhandle.write("\n")

