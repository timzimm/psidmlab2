import h5py
import json
import sys

file = h5py.File(sys.argv[1], 'r')
param_string = file['/'].attrs['parameters'][0].decode("ascii")
p = json.loads(param_string)
print(json.dumps(p, indent=4))
