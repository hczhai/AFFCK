#!/usr/bin/python

# XYZ Format Different Structure Selector
# Code By Huanchen Zhai, Jan. 7, 2015

import sys, getopt, os, re
from math import *

dmax = 0.45
drel = 0.03
deng = 1e9

dmaxr = 0
drelr = 0

def distance(x, y):
    return sqrt(sum([(a - b) ** 2 for a, b in zip(x, y)]))

def im_dist(list):
    k = []
    for i in range(0, len(list)):
        for j in range(i + 1, len(list)):
            k.append(distance(list[i], list[j]))
    k.sort()
    return k

def energy_check(lle, le):
    global deng
    for i in range(0, len(lle)):
        if abs(lle[i] - le) <= deng:
            return 0
    return -1

def similar_check(lld, ld):
    global dmax, drel, dmaxr, drelr
    x = True
    dmaxr = 0
    drelr = 0
    for i in range(0, len(lld)):
        x = False
        y = 0.0
        z = 0.0
        for j in range(0, len(ld)):
            y += abs(lld[i][j] - ld[j])
            z += lld[i][j]
            if abs(lld[i][j] - ld[j]) > dmaxr: dmaxr = abs(lld[i][j] - ld[j])
            if abs(lld[i][j] - ld[j]) > dmax:
                x = True
                break
        if y / z > drelr: drelr = y / z
        if y / z > drel or x:
            x = True
        else:
            return i
    return -1

def read_xyz(xyz):
    xyzl = xyz.split("\n")
    xyzl = [z for z in xyzl if z != ""]
    la = []
    for j in range(2, len(xyzl)):
        al = xyzl[j].split(" ")
        al = [a for a in al if a != ""]
        al = [float(a) for a in al[1:]]
        la.append(al)
    return im_dist(la)

def xyz_comp(lxyz, leng):
    lld, inds, indsn, igs, lle = [], [], [], [], []
    for i in range(0, len(lxyz)):
        fl = open(lxyz[i], "r")
        x = fl.read()
        fl.close()
        ld = read_xyz(x)
        isc = energy_check(lle, leng[i])
        if isc != -1:
            isc = similar_check(lld, ld)
        if isc == -1:
            lld.append(ld)
            lle.append(leng[i])
            inds.append(i)
            indsn.append(1)
        else:
            indsn[isc] += 1
        igs.append(isc)
    return (inds, indsn, igs)

def main_comp(x, y):
    fl = open(x, "r")
    xx = fl.read()
    fl.close()
    fl = open(y, "r")
    yy = fl.read()
    fl.close()
    ldx = read_xyz(xx)
    ldy = read_xyz(yy)
    lld = [ldx]
    isc = similar_check(lld, ldy)
    if isc == -1:
        return 'Different!' 
    else:
        return 'Similar!'

def ds(x, n):
    return "0" * (n - len(str(x))) + str(x)

def main(eng_list_f, inds_col, inp_format, op_f, op_d, ip_u, op_u, op_f2, ren, ap, selt):
    elf = open(eng_list_f, "r")
    el = elf.read().replace("\t", " ").split("\n")
    ell = [[gg for gg in g.split(" ") if gg != ""] for g in el if g != ""]
    elf.close()
    inds = []
    engs = []
    ut = 1.0
    utt = 1.0
    if "h" in ip_u:
        ut = 27.21138505
    if "h" in op_u:
        utt = 1.0 / 27.21138505
    for i in ell:
        m = re.search(r'[0-9]+', i[inds_col[0]])
        inds.append(int(m.group()))
        engs.append(float(i[inds_col[1]]) * ut)
    lxyz = []
    lif = len(inp_format)
    inp_format = inp_format.replace("#", "%d");
    nlif = len(inp_format) - lif
    for i in inds:
        lxyz.append(inp_format % ((i, ) * nlif))
    indsxa, indsxb, igs = xyz_comp(lxyz, engs)
    indsx = (indsxa, indsxb)
    fop = open(op_f, "w")
    if op_f2 != "":
        fop2 = open(op_f2, "w")
        for i in range(0, min(len(igs), selt)):
            fop2.write("%d %d %d (.xyz) %.10f (%.4f) (*%s)\n" % (i, 
                indsxa[igs[i]] if igs[i] != -1 else -1, inds[i], engs[i] * utt, 
                    engs[i] * utt - (engs[0] * utt if ren == None else ren), op_u))
        fop2.close()
    if not os.path.exists(op_d): os.mkdir(op_d)
    if len(indsx[0]) != 0: ind0 = indsx[0][0]
    print "%d / %d" % (len(indsxa), len(lxyz), )
    k = 0
    ffp = open("%s/all.xyz" % op_d, "w")
    for ind, indn in zip(*indsx):
        fop.write("%d (*%d/%d (%.1f%%)) %d (.xyz) %.10f (%.4f) (*%s)\n" % (ind, 
            indn, len(lxyz), indn * 100.0 / len(lxyz), inds[ind], engs[ind] * utt, 
            engs[ind] * utt - (engs[ind0] * utt if ren == None else ren), op_u))
        fn = os.path.basename(lxyz[ind])
        if ap:
            fft = "%s-%s, E = %.10f" % (ds(ind, 4), fn, engs[ind] * utt, )
        else:
            fft = "%s, E = %.10f" % (fn, engs[ind] * utt, )
        fl = open(lxyz[ind], "r")
        x = fl.readlines()
        x = [x[0], fft + "\n"] + x[2:]
        fl.close()
        ffp.write(''.join(x))
        k = k + 1
        if k == selt: break
    ffp.flush()
    ffp.close()
    fop.flush()
    fop.close()


def usage():
    st = """
    XYZ Format Different Structure Selector

    Usage: %s 
                [-h] -l <path> -f <path> [-i <num>] [-e <num>] [-r <num>]
                [-o <path>] [-d <path>] [--dmax=<num>] [--drel=<num>] [--deng=<num>]
                [--input-unit=ev|hartree] [--output-unit=ev|hartree]
                [--prefix=on|off] [--select=<num>]

           %s -c [--dmax=<num>] [--drel=<num>]
                <xyz_file_name_1> <xyz_file_name_2>
                

    >>>> (Part I) Structure Selection
    Options:
        -h:         show help message.
        -l path:    energy list file
        -i num:     index column number (of '-l' file) (default 0)
        -e num:     energy column number (of '-l' file) (default 1)
        -r num:     reference energy in ev (default auto)
        -f path:    input xyz file path (use # to denote file index, 
                    # will be replaced by number got from '-i' column)
        -o path:    output txt file name (default struct_sel.txt)
        -d path:    output dir name (default output)

        --input-unit=ev|hartree:    (default hartree)
        --output-unit=ev|hartree:   (default hartree)
        --output-similar=path:      additional output txt file name
                                    (default None)
        --dmax=num:     max individual diff of dist in similar struct
                        in the input length unit from xyz (default 0.45)
        --drel=num:     max relative diff of dist in similar struct
                        (default 0.03)
        --deng=num:     max energy diff in similar struct in ev
                        (default 1e9, i.e. not applied)
        --prefix=on|off:add prefix ####- to the output xyz files
                        (default on)
        --select=num:   select fisrt <num> structs to output
                        (default 99999, i.e. not applied)

        (Only '-l' and '-f' options are required.)

    Example: %s -l "results/singlet/sorted-singlet.txt"
                     -f "results/singlet/pos_#/posf_#.xyz"
                     --output-unit=ev
                     --output-similar="struct_sim.txt"

    Output file format:
        '-o' file: meanings of columns:
            index | (*occurred times %%) | original index (.xyz) | 
            | absolute energy | (relative energy) (*unit)
        (In '-o' file, each line denotes a unique struct shape.)

        '--output-similar' file: meanings of columns:
            index | (similar with index) | original index (.xyz) |
            | absolute energy | (relative energy) (*unit)
        ('--output-similar' file contains all structs.)

    >>>> (Part II) Structure Comparison
    Options:
        -c:         required to show the task is compare two
                    structures
        xyz_file_name_1:       path to the first xyz file
        xyz_file_name_2:       path to the second xyz file
        --dmax=num:     max individual diff of dist in similar struct
                        in the input length unit from xyz (default 0.45)
        --drel=num:     max relative diff of dist in similar struct
                        (default 0.03)
        --deng=num:     not applied, since there are no energy info provided

    Output will be put on the screen:
        'Similar!' if the two are similar.
        'Different!' if the two are different.

    Example: %s -c pos_0.xyz pos_1.xyz

    """ % (sys.argv[0],sys.argv[0],sys.argv[0],sys.argv[0], ) 
    return st

if __name__ == '__main__':
    eng_list_f, inp_format = [None] * 2
    op_f = "struct_sel.txt"
    op_f2 = ""
    op_d = "output"
    inds_col = [0, 1]
    ip_u = "hartree"
    op_u = "hartree"
    opts, args = getopt.getopt(sys.argv[1:], "hcl:i:e:f:o:d:r:", ["input-unit=", 
        "output-unit=", "output-similar=", "dmax=", "drel=", "deng=", "prefix=", 
        "select="])
    if len(opts) == 0: opts.append(("-h", None))
    comp = False
    ren = None
    ap = True
    selt = 99999
    for op, val in opts:
        if op == "-l":
            eng_list_f = val
        elif op == "-i":
            inds_col[0] = int(val)
        elif op == "-e":
            inds_col[1] = int(val)
        elif op == "-f":
            inp_format = val
        elif op == "-o":
            op_f = val
        elif op == "-d":
            op_d = val
        elif op == "--input-unit":
            ip_u = val
        elif op == "--output-unit":
            op_u = val
        elif op == "--output-similar":
            op_f2 = val
        elif op == "--dmax":
            dmax = float(val)
        elif op == "--drel":
            drel = float(val)
        elif op == "--deng":
            deng = float(val)
        elif op == "--prefix":
            ap = val == "on"
        elif op == "--select":
            selt = int(val)
        elif op == "-c":
            comp = True
        elif op == "-r":
            ren = float(val)
        else:
            print (usage())
            sys.exit()
    if (eng_list_f == None or inp_format == None) and not comp:
        print "Require -l and -f parameters!"
        sys.exit()
    elif comp and len(args) != 2:
        print "Require two file paths!"
        sys.exit()

    if not comp:
        main(eng_list_f, inds_col, inp_format, op_f, op_d, ip_u, op_u, op_f2, ren, ap, selt)
        print "done!"
    else:
        print main_comp(args[0], args[1])
        print "dmax = %.4f, drel = %.4f" % (dmaxr, drelr)

    sys.exit()
