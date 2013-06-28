import numpy as np
import os

class group_tab:
  def __init__(self, basedir, snapnum):
    self.filebase = basedir + "/groups_" + str(snapnum).zfill(3) + "/group_tab_" + str(snapnum).zfill(3) + "."

    done_flag=False
    filenum=0
    while (done_flag==False):	
	
	f=open(self.filebase+str(filenum), "rb")
	Ngroups=np.fromfile(f, count=1, dtype="int32")[0]
	self.TotNgroups=np.fromfile(f, count=1, dtype="int32")[0]
	Nids=np.fromfile(f, count=1, dtype="int32")[0]
	self.TotNids=np.fromfile(f, count=1, dtype="uint64")[0]
	self.NTask=np.fromfile(f, count=1, dtype="int32")[0]

	GroupLen=np.fromfile(f, count=Ngroups, dtype="int32")
	GroupOffset=np.fromfile(f, count=Ngroups, dtype="int32")
	GroupMass=np.fromfile(f, count=Ngroups, dtype="float32")
	GroupCM=np.fromfile(f, count=3*Ngroups, dtype="float32").reshape(Ngroups,3)
	GroupVel=np.fromfile(f, count=3*Ngroups, dtype="float32").reshape(Ngroups,3)
	GroupLenType=np.fromfile(f, count=6*Ngroups, dtype="int32").reshape(Ngroups,6)
	GroupMassType=np.fromfile(f, count=6*Ngroups, dtype="float32").reshape(Ngroups,6)
	GroupSfr=np.fromfile(f, count=Ngroups, dtype="float32")

        if (filenum==0):
                self.GroupLen=GroupLen
                self.GroupOffset=GroupOffset
                self.GroupMass=GroupMass
                self.GroupCM=GroupCM
                self.GroupVel=GroupVel
                self.GroupLenType=GroupLenType
                self.GroupMassType=GroupMassType
                self.GroupSfr=GroupSfr
        else:
                self.GroupLen=np.append(self.GroupLen, GroupLen)
		self.GroupOffset=np.append(self.GroupOffset, GroupOffset)
		self.GroupMass=np.append(self.GroupMass, GroupMass)
		self.GroupCM=np.append(self.GroupCM, GroupCM, axis=0)
		self.GroupVel=np.append(self.GroupVel, GroupVel, axis=0)
		self.GroupLenType=np.append(self.GroupLenType, GroupLenType, axis=0)
		self.GroupMassType=np.append(self.GroupMassType, GroupMassType, axis=0)
		self.GroupSfr=np.append(self.GroupSfr, GroupSfr)

	curpos = f.tell()
        f.seek(0,os.SEEK_END)
        if curpos != f.tell(): print "Warning: finished reading before EOF for file",filenum


	f.close()
	
	filenum+=1

	if (filenum==self.NTask):
		done_flag=True


class group_ids:
  def __init__(self, basedir, snapnum, longids=False):
    self.filebase = basedir + "/groups_" + str(snapnum).zfill(3) + "/group_ids_" + str(snapnum).zfill(3) + "."

    done_flag=False
    filenum=0
    while (done_flag==False):

        f=open(self.filebase+str(filenum), "rb")
        Ngroups=np.fromfile(f, count=1, dtype="int32")[0]
        self.TotNgroups=np.fromfile(f, count=1, dtype="int32")[0]
        Nids=np.fromfile(f, count=1, dtype="int32")[0]
        self.TotNids=np.fromfile(f, count=1, dtype="uint64")[0]
        self.NTask=np.fromfile(f, count=1, dtype="int32")[0]
	Nprevids=np.fromfile(f, count=1, dtype="int32")[0]

	if longids:
		IDs=np.fromfile(f, count=Nids, dtype="uint64")
	else:
		IDs=np.fromfile(f, count=Nids, dtype="uint32")

        if (filenum==0):
		self.IDs=IDs
        else:
		self.IDs=np.append(self.IDs, IDs)
        f.close()

        filenum+=1

        if (filenum==self.NTask):
                done_flag=True


