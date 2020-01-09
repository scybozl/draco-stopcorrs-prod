import os
import sys

path = sys.argv[1]
cwd = os.getcwd()
l = len(os.listdir(path))

nfail = 0

faillist = []
lists = []

for i,dirs in enumerate(os.listdir(path)):

  curdir = cwd+'/'+path+dirs
  os.chdir(curdir)
  hand = dirs.split('-')[0]
  PRmss = float(dirs.split('-Mst.')[1].split('-M')[0])
  PRcim = float(dirs.split('-M1.')[1].split('-mt')[0])

  print i+1, "/", l, " ---> Cross-sections for ", hand, ", mstop = ", PRmss, ", mchi = ", PRcim
  print " -----------------------------------------------------------------------------"

  if not os.path.exists("xsec.out"):
    print " No xsec file!"
    print " #############################################################################\n"
    continue

  xsecs = open("xsec.out", 'r')
  lines = xsecs.readlines()
  xsecemu0 = lines[1].split('\n')[0]
  xsecemu1 = lines[4].split('\n')[0]
  xsecmue0 = lines[7].split('\n')[0]
  xsecmue1 = lines[10].split('\n')[0]

  print "    e+, mu- (+1j) = ", xsecemu0, " ( ", xsecemu1, " ) ", "  \n -- e-, mu+ (+1j) = ", xsecmue0, " ( ", xsecmue1, " ) "

  if float(xsecemu0) > 1e-5 and float(xsecemu1) > 1e-5 and float(xsecmue0) > 1e-5 and float(xsecmue1) > 1e-5:
    print "    -> \033[92m OK \033[0m"
    lists += [[hand, PRmss, PRcim, 1]]
  else:
    nfail += 1
    print "    -> \033[91m FAIL \033[0m"
    lists += [[hand, PRmss, PRcim, 0]]
    faillist += [[hand, PRmss, PRcim, 0]]

  print " #############################################################################\n"
  xsecs.close()

realfail = 0
realfaillist = []
for e in lists:

  if e[3] == 1: continue
  if e[0] == 'L': handto = 'R'
  if e[0] == 'R': handto = 'L'

  for j in lists:
    if e[1] == j[1] and e[2] == j[2] and j[0] == handto:
      if j[3] == 0: 
        realfail += 1
        realfaillist += [str(j[1]) + " / " + str(j[2])]

print "\n\n Number of divergent ( i.e. tiny ) PS points: ", nfail, " / ", l
print " Number of double-handed divergent ( i.e. tiny ) PS points: ", realfail, " / ", l

print " ------ Bad points ------"
for i in list(set(realfaillist)):
  print "   ", i
