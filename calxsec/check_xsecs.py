import os
import sys

path = sys.argv[1]
cwd = os.getcwd()
l = len(os.listdir(path))

nfail = 0

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
  else:
    nfail += 1
    print "    -> \033[91m FAIL \033[0m"

  print " #############################################################################\n"
  xsecs.close()

print "\n\n Number of divergent ( i.e. tiny ) PS points: ", nfail, " / ", l

