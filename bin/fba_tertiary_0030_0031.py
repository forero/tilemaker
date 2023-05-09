# generic python
import tempfile
import shutil
import os
import numpy as np
import argparse

# from the astro comunity
import healpy

# from astropy
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy import units as u

# from desihub
from desitarget.targets import encode_targetid
from fiberassign.fba_tertiary_io import assert_tertiary_targ
from fiberassign.fba_tertiary_io import assert_tertiary_prio
from desisurvey.tileqa import lb2uv
from desiutil import dust

def create_tertiary_priority_dict():
    return {
        "PAL5_GAIA" : {"PRIORITY_INIT" : 8000, "NGOAL" : 4},
        "PAL5_NOGAIA" : {"PRIORITY_INIT" : 7000, "NGOAL" : 4},
        "SGR_GAIA" : {"PRIORITY_INIT" : 6000, "NGOAL" : 4},
        "SGR_NOGAIA" : {"PRIORITY_INIT" : 5000, "NGOAL" : 4},
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
    d.write(outfn, overwrite=True)

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

def prepare_data_table(d, prognum, obsconds, goaltime, np_rand_seed):
    d["TARGETID"] = encode_targetid(release=8888, brickid=prognum, objid=np.arange(len(d)))
    d["SUBPRIORITY"] = np.random.uniform(size=len(d))
    d.meta["EXTNAME"] = "TARGETS"
    d.meta["FAPRGRM"] = "tertiary{}".format(prognum)
    d.meta["OBSCONDS"] = obsconds
    d.meta["SBPROF"] = "PSF"
    d.meta["GOALTIME"] = goaltime
    d.meta["RANDSEED"] = np_rand_seed
    return d

def create_targets_file(output_path, prognum, obsconds, goaltime, input_filename):
    outfn = os.path.join(output_path, "tertiary-targets-{:04d}.fits".format(prognum))
    np_rand_seed = 1234
    np.random.seed(np_rand_seed)
    
    d = read_input_data(input_filename)
    d.meta["INPUTFN"] = input_filename
    d = prepare_data_table(d, prognum, obsconds, goaltime, np_rand_seed)
    d.write(outfn, overwrite=True)


# AR grid dither of ~1 patrol radius
# AR fiber patrol diameter = 12mm
# AR                       = 0.047 deg with taking 70.4 um/arcsec
# AR we take 0.048 deg
def get_centers(field_ra, field_dec, npt=5, rad=0.048):
    vals = rad * np.arange(-10, 11)
    if npt > vals.size ** 2:
        raise ValueError("get_centers: npt={} is larger than {} -> need to re-write get_centers!".format(npt, vals.size ** 2))
    ras, decs = [], []
    for i in range(len(vals)):
        ras += (field_ra + vals[i] / np.cos(np.radians(field_dec + vals))).tolist()
        decs += (field_dec + vals).tolist()
    ras, decs = np.array(ras), np.array(decs)
    field_c = SkyCoord(field_ra * u.degree, field_dec * u.degree, frame="icrs")
    cs = SkyCoord(ras * u.degree, decs * u.degree, frame="icrs")
    ii = cs.separation(field_c).value.argsort()
    ras, decs = ras[ii], decs[ii]
    ras, decs = ras[:npt], decs[:npt]
    ras[ras < 0] += 360
    return ras, decs

# AR adapted from https://github.com/desihub/desisurvey/blob/94b02bdae04137526bf98dcec0dca8bd29a231d3/py/desisurvey/tileqa.py#L525-L554 
def get_ebv_meds(tiles):
    nside = 512
    theta, phi = healpy.pix2ang(nside, np.arange(12*nside**2))
    la, ba = phi*180./np.pi, 90-theta*180./np.pi
    from desiutil import dust
    ebva = dust.ebv(la, ba, frame='galactic',
                    mapdir=os.getenv('DUST_DIR')+'/maps', scaling=1)
    if isinstance(tiles, Table):
        ra = tiles['RA'].data
        dec = tiles['DEC'].data
    else:
        ra = tiles['RA']
        dec = tiles['DEC']
    coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg,
                     frame='icrs')
    coordgal = coord.galactic
    lt, bt = coordgal.l.value, coordgal.b.value
    uvt = lb2uv(lt, bt)
    fprad = 1.605
    ebv_meds = np.zeros(len(tiles))
    for i in range(len(tiles)):
        ind = healpy.query_disc(nside, uvt[i], fprad*np.pi/180.)
        ebv_meds[i] = np.median(ebva[ind])
    return ebv_meds

def create_tiles_table(tileids, ras, decs, program):
    d = Table()
    d["TILEID"] = tileids
    d["RA"], d["DEC"] = ras, decs
    d["PROGRAM"] = program.upper()
    d["IN_DESI"] = True
    d["EBV_MED"] = get_ebv_meds(d)
    d["DESIGNHA"] = 0  # eddie s email from 9/6/2022 10:54am
    return d

def create_tiles_file(output_path, prognum, program, tileids):
    ntile = len(tileids)
    field_ra, field_dec = 229.5, 1.0

    ras, decs = get_centers(field_ra, field_dec, npt=ntile)
    ras, decs = ras.round(3), decs.round(3)

    d = create_tiles_table(tileids, ras, decs, program)
    
    outfn = os.path.join(output_path, "tertiary-tiles-{:04d}.ecsv".format(prognum))
    d.write(outfn, overwrite=True)

def assert_files(output_path, prognum):
    targetsfn = os.path.join(output_path, "tertiary-targets-{:04d}.fits".format(prognum))
    prioritiesfn = os.path.join(output_path, "tertiary-priorities-{:04d}.ecsv".format(prognum))

    assert_tertiary_prio(prognum, prioritiesfn, targetsfn)
    assert_tertiary_targ(prognum, targetsfn)


def main(out_dir="./"):
    
    
    # Step 1.
    # Copy the input target file into the working directory directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    src_file_path = "/global/cfs/cdirs/desi/users/raichoor/fiberassign-desi2-pal5/v0/inputcats/Pal5andSgr_desi2sel.fits"

    # Construct the destination file path inside the output directory
    dst_file_path = os.path.join(out_dir, os.path.basename(src_file_path))

    # Copy the file from location A to the output directory
    shutil.copy(src_file_path, dst_file_path)

    print(f"Output directory: {out_dir}")
    print(f"File copied to: {dst_file_path}")
    
    # Step 2. Create the priorities.
    create_priorities_file(out_dir, 30)
    create_priorities_file(out_dir, 31)

    # Step 3. Create the targets
    create_targets_file(out_dir, 30, "BRIGHT", 600, dst_file_path)
    create_targets_file(out_dir, 31, "DARK", 1800, dst_file_path)

    # Step 4. Create the tile files
    create_tiles_file(out_dir, 30, "BRIGHT", [83355, 83356, 83357, 83358, 83359])
    create_tiles_file(out_dir, 31, "DARK", [83360, 83361, 83362, 83363, 83364])

    # Step 5. Assert priorities and target files
    assert_files(out_dir, 30)
    assert_files(out_dir, 31)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create priorities, tiles and targets for fiberassign from an input target file.')
    parser.add_argument('out_dir', type=str, help='Output path to write the priorities, tiles and targets')
    
    args = parser.parse_args()
    main(out_dir=args.out_dir)
