#! /usr/bin/env python
# 
#    beamerstrip.py: a python function to strip a LaTeX beamer class file
#    Outputs a article class .tex file
#    Usage: $ beamerstrip.py beamer.tex article.tex
#    Copyright (C) 2020  Mark E. Fuller
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#    Mark E. Fuller: mark.e.fuller@gmx.de


#setup terminal output later:
#    beamerstrip.py  Copyright (C) 2020  Mark E. Fuller
#    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#    This is free software, and you are welcome to redistribute it
#    under certain conditions; type `show c' for details.


################################################################################


# Import stuff.
import sys

beamerfile = sys.argv[1]#"slides.tex"
articlefile = sys.argv[2]#"mdtest.tex"

#read in beamer file (bf)
with open(beamerfile, "r") as bf:
    bflines = bf.readlines()

#init output (article) file (af)
aflines=[]

#first write the preamble information that we want here
# use raw strings (r"...") to avoid problems with backslash
aflines.append(r"\documentclass[a4paper,10pt]{article}")
aflines.append("")
aflines.append(r"\usepackage[utf8]{inputenc}")
aflines.append(r"%\usepackage[ngerman]{babel}        % Deutsch, neue Rechtschreibung")
aflines.append(r"\usepackage[english]{babel}       % English")
aflines.append("")
aflines.append(r"\usepackage[T1]{fontenc}           % Font encoding (don't change!)")
aflines.append(r"\usepackage{lmodern}               % Select Linux Modern Fonts (don't change)")
aflines.append(r"\usepackage{sansmathfonts}         % Sans fonts in math environments")
aflines.append(r"\usepackage{textcomp}              % fix 'missing font symbols' warning")
aflines.append(r"\renewcommand{\rmdefault}{phv}     % Arial like (Helvetica)")
aflines.append(r"\renewcommand{\sfdefault}{phv}     % Arial like (Helvetica)")
aflines.append("")
aflines.append(r"\usepackage{graphicx}              % needed to include graphics (don't change)")
aflines.append(r"\usepackage{epstopdf}              % required to include eps files")
aflines.append(r"\usepackage[encoding,filenameencoding=utf8]{grffile} % allow utf8 file names in graphics")
aflines.append(r"\usepackage{listings}                           % for lstlisting and \lstinline|..|")
aflines.append(r"\usepackage{hyperref}")
aflines.append(r"\usepackage{amsmath}")
aflines.append(r"\usepackage{amssymb}")
aflines.append(r"\usepackage{sansmath}")
aflines.append(r"\usepackage{tabularx}")
aflines.append(r"\usepackage{booktabs}")
aflines.append(r"\usepackage[version=3]{mhchem} % Formula subscripts using \ce{}")


#possible additional things:
#flag \cite, \footcite, \footfullcite... for replacement with link ?
#flag images with replacement with wiki location of file ?

AutoAppend = False # do not append lines by default
for line in bflines:
#use strip() to remove whitespace
    if line.strip().startswith(r"%"): #comment
        continue # skip writing in the comments for now
    elif line.strip().startswith(r"\title"):
        aflines.append(line.strip())
    elif line.strip().startswith(r"\author"):
        aflines.append(line.strip())
    elif line.strip().startswith(r"\begin{document}"):
        aflines.append(line.strip())
        aflines.append(r"\maketitle") #not in beamer file
        AutoAppend = True #this is the start of stuff to copy
    elif line.strip().startswith(r"\bib"): #end of file bibliography
        AutoAppend = False
        aflines.append(r"\end{document}")
        break # no need to go any further
    elif line.strip().startswith(r"\end{document}"):
        aflines.append(line.strip())
        break
    elif line.strip().startswith(r"\begin{frame}"):
        aflines.append(line.strip().replace(r"\begin{frame}",r"\subsubsection"))
    elif line.strip().startswith(r"\end{frame}"):
        aflines.append("") #replace with break after frame
    elif line.strip()=="": #blank line
        continue #skip this line: I'm already reformatting, so I don't need original blank lines
    else:
        if AutoAppend:
            aflines.append(line.strip())

#write update file    
with open(articlefile, "w") as af:
    #write the single string creating by joining all the strings in the aflines list
    #join takes all the strings in list and merges to one string, separating with newline('\n')
    af.write('\n'.join(aflines))
