#HDF5 library wrapper for:
#*PyTables (http://www.pytables.org)
#*h5py (http://code.google.com/p/h5py/)
#
# Mark Vogelsberger (mvogelsb@cfa.harvard.edu)
import sys

try:
	import tables
	use_tables=True
except ImportError:
	import h5py
	use_tables=False


def OpenFile(fname, mode = "r"):
	if (use_tables):
		return tables.openFile(fname, mode = mode)
	else:
		return h5py.File(fname, mode)

def GetData(f, dname):
	if (use_tables):
		return f.root._f_getChild(dname)
	else:
		return f[dname]

def GetGroup(f, gname):
	if (use_tables):
		return f.root._f_getChild(gname) 
	else:
		return f[gname]

def GetAttr(f, gname, aname):
	if (use_tables):
		return f.root._f_getChild(gname)._f_getAttr(aname) 
	else:
		 return f[gname].attrs[aname]

		
def Contains(f, gname, cname):
	if (use_tables):
		if gname=="":
			return f.root.__contains__(cname) 
		else:
			return f.root._f_getChild(gname).__contains__(cname)
	else:
		if gname=="":
			return f.__contains__(cname)
		else:
			return f[gname].__contains__(cname)


def ContainsGroup(group, cname):
	return group.__contains__(cname)

		
def CreateArray(f, where, aname, aval):
	if (use_tables):
		f.createArray(where, aname, aval)
	else:
		print "NOT IMPLEMENTED!"
		sys.exit()


def CreateGroup(f, gname):
	if (use_tables):
		return f.createGroup(f.root, gname)
	else:
		print "NOT IMPLEMENTED!"
		sys.exit()


def SetAttr(where, aname, aval):
	if (use_tables):
		setattr(where._v_attrs, aname, aval)
	else:
		print "NOT IMPLEMENTED!"
		sys.exit()
