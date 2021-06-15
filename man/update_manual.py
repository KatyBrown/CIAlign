#!/usr/bin/env python
import os

'''
Adds the default, minimum and maximum values for the cleaning functions
to the user guide template.
Converts to pdf
Adds the correct header to the GitHub readme
'''

def main():
    ci_dir = "../CIAlign"
    ranges = [line.strip().split("\t")
              for line in open("%s/ranges.txt" % ci_dir)]
    # Defaults
    defs = {x[0]: x[1] for x in ranges}
    # Minima
    minis = {x[0]: x[2] for x in ranges}
    # Maxima
    maxis = {x[0]: x[3] for x in ranges}
    
    statement = 'cat user_guide_template.md'

    for de in defs:
        statement += " | sed 's/%s_def/%s/g'" % (de, defs[de])
        statement += " | sed 's/%s_min/%s/g'" % (de, minis[de].replace("/", "\/"))
        statement += " | sed 's/%s_max/%s/g'" % (de, maxis[de].replace("/", "\/"))
    
    statement += "> user_guide.md"
    os.system(statement)
    
    statement = """\
    pandoc --from markdown --to latex -o user_guide.pdf user_guide.md \
        --template=template.tex
    echo '# CIAlign' > ../README.md
    awk 'NR > 7' user_guide.md >> ../README.md"""
    os.system(statement)


if __name__ == "__main__":
    main()