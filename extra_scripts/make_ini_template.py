#!/usr/bin/env python3
"""
Converts the output of CIAlign --help into a config file with the
default options
"""
import sys

def main():

    infile = sys.argv[1] # python3 CIAlign --help > temp.txt
    outfile = sys.argv[2] # ini file

    pgroups = [line.strip().split("\t")
               for line in open("extra_scripts/param_groups.txt").readlines()]
    groupD = dict()
    for group, params in pgroups:
        params = params.split(",")
        for param in params:
            groupD[param] = group
    lines = [line.strip() for line in open(infile).readlines()]
    bigstring = " ".join(lines)
    bigstring = bigstring.split("Required Arguments: ")[1]
    bigstring = bigstring.split("-h")[0]
    arguments = bigstring.split("--")[1:-1]
    out = open(outfile, "w")
    out.write("# CIAlign Parameters\n")
    groups_done = set()
    for argument in arguments:
        nam = argument.split(" ")[0]
        if nam in groupD:
            group = groupD[nam]
            if group not in groups_done:
                out.write("\n[%s]\n" % group)
                groups_done.add(group)
            parts = argument.split(" ")[1:]
            if parts[0].upper() == parts[0]:
                parts = parts[1:]
            desc = " ".join(parts)
            desc = desc.replace("  ", " ").strip()
            desc = desc.replace("Optional Arguments:", "").replace(" -h,", "")
            if "Required" not in desc:
                if "Default: " not in desc:
                    print (desc)
                desc, d2 = desc.split("Default: ")
            else:
                desc = desc.replace(" Required", "")
                d2 = "example.py"
            if "config" not in desc:
                if d2 != "None":
                    out.write("# %s\n%s = %s\n" % (desc, nam, d2))
                else:
                    out.write("# %s\n%s =\n" % (desc, nam))
    out.close()

if __name__ == "__main__":
    main()
