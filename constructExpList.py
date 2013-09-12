import numpy as np

noutputs = 1024
expfacti = 0.0212765957447

snapshots = range(0,noutputs)
expfact = np.linspace(expfacti,1,noutputs)

f = open("ExpansionList_" + str(noutputs),'w')

for i in xrange(0,len(snapshots)):
    outexpt = '%0.9f' % (expfact[i],)
    f.write(str(outexpt) + " 1\n")

f.close()