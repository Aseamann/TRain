import chimera
import os
import Midas
import math

directory = os.getcwd()
# path = directory + "6m0j.pdb"
energy_table = "test.csv"
chimera.runCommand("split #0 chains")


def floatRgb(mag, cmin, cmax):
    """ Return a tuple of floats between 0 and 1 for R, G, and B. """
    # Normalize to 0-1
    try:
        x = float(mag-cmin)/(cmax-cmin)
    except ZeroDivisionError:
        x = 0.5  # cmax == cmin
    blue = min((max((4*(0.75-x), 0.)), 1.))
    red = min((max((4*(x-0.25), 0.)), 1.))
    green = min((max((4*math.fabs(x-0.5)-1., 0.)), 1.))
    return red, green, blue


def rgb(mag, cmin, cmax):
    """ Return a tuple of integers, as used in AWT/Java plots. """
    red, green, blue = floatRgb(mag, cmin, cmax)
    return int(red*255), int(green*255), int(blue*255)


def strRgb(mag, cmin, cmax):
    """ Return a hex string, as used in Tk plots. """
    return "#%02x%02x%02x" % rgb(mag, cmin, cmax)


def get_mag(current, min_i, max_i):
    return (1.0 / (max_i - min_i)) * abs(current)


with open(energy_table, "r") as file:
    values = []
    ran = {}
    min_i = 10000000.0
    max_i = -10000000.0
    for line in file:
        split = line.split(",")
        if split[0] != "ACE2":
            # 0: ACE2 AA, 1: value, 2: Spike AA
            values.append([split[0], float(split[1]), split[2][:-1]])
    for each in values:
        if each[1] > max_i:
            max_i = each[1]
        if each[1] < min_i:
            min_i = each[1]
    for each in values:
        if each[1] <= 0.0:
            if each[0] in ran.keys():
                if ran[each[0]] > each[1]:
                    ran[each[0]] = each[1]
                    # Colors ACE2
                    chimera.runCommand(
                        "color " + strRgb(get_mag(each[1], min_i, 0), 0, 1) + ",r,s #0.1:" + str(each[0].split(" ")[1]))
                    # Colors Spike
                    chimera.runCommand(
                        "color " + strRgb(get_mag(each[1], min_i, 0), 0, 1) + ",r,s #0.2:" + str(each[2].split(" ")[1]))
            else:
                ran[each[0]] = each[1]
                # Colors ACE2
                chimera.runCommand("color " + strRgb(get_mag(each[1], min_i, 0), 0, 1) + ",r,s #0.1:"
                                   + str(each[0].split(" ")[1]))
                # Shows the side chain of ACE2
                chimera.runCommand("display #0.1:" + str(each[0].split(" ")[1]))
                # Colors Spike
                chimera.runCommand(
                        "color " + strRgb(get_mag(each[1], min_i, 0), 0, 1) + ",r,s #0.2:" + str(each[2].split(" ")[1]))
                # Shows the side chain of spike
                chimera.runCommand("display #0.2:" + str(each[2].split(" ")[1]))

