import os
import pylab
import lsst.sims.atmosphere.transmission.modtranCards as mc
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass

# Instantiate object
m = mc.ModtranCards()

# Read default card template
m.readCardTemplate()
#print m._cardTemplate

# Read default parameter formats
m.readParameterFormats()
#print m._paramFormats

# Read a test input parameter file (in Michel's format)
paramList = m.readParamValues_M('tmp2_parmlist.dat')
print paramList
#for run in paramList:
#    print run['ICLD']

# Write modtran input cards
m.writeModtranCards(paramList[0], 'tmp2')

# Run modtran on this file.
m.runModtran('tmp2')

atm = Bandpass()
atm.readThroughput('tmp2.psc')
pylab.plot(atm.wavelen, atm.sb, 'b-')
pylab.show()

m.cleanModtran('tmp2')

exit()

# compare tmp.tp5 with tmp.tp5_ref
file1 = open('tmp.tp5', 'r')
file2 = open('ref_tmp.tp5', 'r')
for line1, line2 in zip(file1, file2):
    print 'new', line1.strip()
    print 'ref', line2.strip()
file1.close()
file2.close()
