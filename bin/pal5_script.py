import tempfile
import shutil
import os
import numpy as np
from astropy.table import Table, vstack
from desitarget.targets import encode_targetid
from fiberassign.fba_tertiary_io import assert_tertiary_targ

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

def create_priorities_file(output_path, prognum):
    outfn = os.path.join(output_path, "tertiary-priorities-{:04d}.ecsv".format(prognum))

    mydict = create_tertiary_priority_dict()
    myd = create_empty_data_dict()
    populate_data_dict(mydict, myd)
    d = create_table_from_dict(myd)
    d.write(outfn)

def read_input_data(filename):
    d = Table.read(filename)
    d.meta = None
    d["TERTIARY_TARGET"] = d["SAMPLE"].astype(str)
    for key in ["ID", "FLUX_G", "FLUX_R"]:
        d["ORIG_{}".format(key)] = d[key]
    for key in ["PMRA", "PMDEC", "REF_EPOCH"]:
        d[key] = d[key].astype(np.float32)
    d.keep_columns(["RA", "DEC", "PMRA", "PMDEC", "REF_EPOCH", "TERTIARY_TARGET", "CHECKER", "ORIG_ID", "ORIG_FLUX_G", "ORIG_FLUX_R"])
    sel = (~np.isfinite(d["PMRA"])) | (~np.isfinite(d["PMDEC"]))
    sel |= (d["PMRA"] == -999.) | (d["PMDEC"] == -999.)
    d["PMRA"][sel], d["PMDEC"][sel], d["REF_EPOCH"][sel] = 0., 0., 2015.5
    return d

def prepare_data_table(d, prognum, np_rand_seed):
    d["TARGETID"] = encode_targetid(release=8888, brickid=prognum, objid=np.arange(len(d)))
    d["SUBPRIORITY"] = np.random.uniform(size=len(d))
    d.meta["EXTNAME"] = "TARGETS"
    d.meta["FAPRGRM"] = "tertiary{}".format(prognum)
    d.meta["OBSCONDS"] = "BRIGHT"
    d.meta["SBPROF"] = "PSF"
    d.meta["GOALTIME"] = 600
    d.meta["RANDSEED"] = np_rand_seed
    return d

def create_targets_file(output_path, prognum, input_filename):
    outfn = os.path.join(output_path, "tertiary-targets-{:04d}.fits".format(prognum))
    np_rand_seed = 1234
    np.random.seed(np_rand_seed)
    
    d = read_input_data(input_filename)
    d.meta["INPUTFN"] = input_filename
    d = prepare_data_table(d, prognum, np_rand_seed)
    d.write(outfn)

    assert_tertiary_targ(prognum, outfn)

def main():
    # Step 1.
    # Copy the input target file into a temporary directory

    # Define the source file path (location A)
    src_file_path = "/global/cfs/cdirs/desi/users/raichoor/fiberassign-desi2-pal5/v0/inputcats/Pal5andSgr_desi2sel.fits"

    # Create a temporary directory
    tmp_dir = tempfile.mkdtemp()

    # Construct the destination file path inside the temporary directory
    dst_file_path = os.path.join(tmp_dir, os.path.basename(src_file_path))

    # Copy the file from location A to the temporary directory
    shutil.copy(src_file_path, dst_file_path)

    print(f"Temporary directory: {tmp_dir}")
    print(f"File copied to: {dst_file_path}")
    
    
    # Step 2. Create the priorities.
    create_priorities_file(tmp_dir, 9998)
    create_priorities_file(tmp_dir, 9999)

    # Step 3. Create the targets
    create_targets_file(tmp_dir, 9998, dst_file_path)
    create_targets_file(tmp_dir, 9999, dst_file_path)

    
    
if __name__ == "__main__":
    main()




# Deleting the temporary file
#shutil.rmtree(tmp_dir)