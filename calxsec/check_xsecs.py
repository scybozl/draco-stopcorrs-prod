import os
import sys

path = sys.argv[1]
cwd = os.getcwd()
l = len(os.listdir(path))

nfail = 0

faillist = []
lists = []

for i,dirs in enumerate(os.listdir(path)):

  if dirs.find('yodas') != -1: continue
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
    badfile = open("status.prec.xs", 'w')
    badfile.write("FAILED")
    badfile.close()
    continue

  xsecs = open("xsec.out", 'r')
  lines = xsecs.readlines()
  xsecemu0 = lines[1].split('(')[0]
  d_xsecemu0 = lines[1].split('+/-')[1].split(')')[0]
  xsecemu1 = lines[4].split('(')[0]
  d_xsecemu1 = lines[4].split('+/-')[1].split(')')[0]
  xsecmue0 = lines[7].split('(')[0]
  d_xsecmue0 = lines[7].split('+/-')[1].split(')')[0]
  xsecmue1 = lines[10].split('(')[0]
  d_xsecmue1 = lines[10].split('+/-')[1].split(')')[0]
#  xsecemu0 = lines[1].split('\n')[0]
#  xsecemu1 = lines[4].split('\n')[0]
#  xsecmue0 = lines[7].split('\n')[0]
#  xsecmue1 = lines[10].split('\n')[0]

  print "    e+, mu- (+1j) = ", xsecemu0, " +/- ", d_xsecemu0, " ( ", xsecemu1, " +/- ", d_xsecemu1, ") ", \
        "  \n -- e-, mu+ (+1j) = ", xsecmue0, " +/- ", d_xsecmue0, " ( ", xsecmue1, " +/- ", d_xsecmue1, ") "

#  If cross-section is unnaturally small
#  if float(xsecemu0) > 1e-5 and float(xsecemu1) > 1e-5 and float(xsecmue0) > 1e-5 and float(xsecmue1) > 1e-5:
  if (float(d_xsecemu0) / float(xsecemu0) < 0.05 and float(d_xsecemu1) / float(xsecemu1) < 0.05 and
      float(d_xsecmue0) / float(xsecmue0) < 0.05 and float(d_xsecmue1) / float(xsecmue1) < 0.05):
    print "    -> \033[92m OK \033[0m"
    lists += [[hand, PRmss, PRcim, 1]]
    okfile = open("status.prec.xs", 'w')
    okfile.write("GOOD")
    okfile.close()
  else:
    nfail += 1
    print "    -> \033[91m FAIL \033[0m"
    lists += [[hand, PRmss, PRcim, 0]]
    faillist += [[hand, PRmss, PRcim, 0]]
    badfile = open("status.prec.xs", 'w')
    badfile.write("FAILED")
    badfile.close()

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

dataL = []
dataR = []
for i in lists:
  if i[0] == 'L':
    if i[3] == 0:
      dataL += [",".join([str(i[1]), str(i[2]), 'red'])]
    elif i[3] == 1:
      dataL += [",".join([str(i[1]), str(i[2]), 'green'])]
  elif i[0] == 'R':
    if i[3] == 0:
      dataR += [",".join([str(i[1]), str(i[2]), 'red'])]
    elif i[3] == 1:
      dataR += [",".join([str(i[1]), str(i[2]), 'green'])]

os.chdir(cwd)

outfile = open('badpoints-L.dat', 'w')
outfile.writelines("\n".join(dataL))
outfile.close()

outfile2 = open('badpoints-R.dat', 'w')
outfile2.writelines("\n".join(dataR))
outfile2.close()
