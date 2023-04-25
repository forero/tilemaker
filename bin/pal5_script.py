import tempfile
import shutil
import os

# Step 1.
# Copy the target file into a temporary directory

# Define the source file path (location A)
src_file_path = "/global/cfs/cdirs/desi/users/raichoor/fiberassign-desi2-pal5/v0/inputcats/Pal5andSgr_desi2sel.fits"

# Create a temporary directory
tmp_dir = tempfile.mkdtemp()

# Construct the destination file path inside the temporary directory
dst_file_path = os.path.join(tmp_dir, os.path.basename(src_file_path))

# Copy the file from location A to the temporary directory
shutil.copy(src_file_path, dst_file_path)

# Do any operations you need with the copied file
# ...

print(f"Temporary directory: {tmp_dir}")
print(f"File copied to: {dst_file_path}")


import numpy as np
from astropy.table import Table, vstack

def create_tertiary_priority_dict():
    return {
        "BRIGHT_PM_BOO3": {"PRIORITY_INIT": 8000, "NGOAL": 4},
        "FAINT_NOPM_BOO3": {"PRIORITY_INIT": 7000, "NGOAL": 4},
        "BRIGHT_PM_GC": {"PRIORITY_INIT": 6000, "NGOAL": 4},
        "FAINT_NOPM_GC": {"PRIORITY_INIT": 5000, "NGOAL": 4},
        "FILLER": {"PRIORITY_INIT": 4000, "NGOAL": 4},
    }

def create_empty_data_dict():
    return {
        key: []
        for key in ["TERTIARY_TARGET", "NUMOBS_DONE_MIN", "NUMOBS_DONE_MAX", "PRIORITY"]
    }

def populate_data_dict(mydict, myd):
    for target in mydict:
        prio_init, ngoal = mydict[target]["PRIORITY_INIT"], mydict[target]["NGOAL"]
        for nmin in range(ngoal + 1):
            if nmin == ngoal:
                nmax, prio = 99, 1
            else:
                nmax, prio = nmin, prio_init + nmin
            myd["TERTIARY_TARGET"].append(target)
            myd["NUMOBS_DONE_MIN"].append(nmin)
            myd["NUMOBS_DONE_MAX"].append(nmax)
            myd["PRIORITY"].append(prio)

def create_table_from_dict(myd):
    d = Table()
    for key in myd:
        d[key] = myd[key]
    return d

def create_priorities_file(path, prognum):
    outfn = os.path.join(path, "tertiary-priorities-{:04d}.ecsv".format(prognum))

    mydict = create_tertiary_priority_dict()
    myd = create_empty_data_dict()
    populate_data_dict(mydict, myd)
    d = create_table_from_dict(myd)
    d.write(outfn)

    
create_priorities_file(tmp_dir, 9998)
create_priorities_file(tmp_dir, 9999)

# Deleting the temporary file
#shutil.rmtree(tmp_dir)