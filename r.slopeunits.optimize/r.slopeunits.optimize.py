#!/usr/bin/env python3
#
############################################################################
#
# MODULE:       r.slopeunits.create for GRASS 8
# AUTHOR(S):    Ivan Marchesini, Massimiliano Alvioli, Carmen Tawalika
# PURPOSE:      Calculate metrics for slope units
# COPYRIGHT:    (C) 2004-2024 by the GRASS Development Team
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
#############################################################################
#
# %module
# % description: Calculate metrics for slope units
# % keywords: raster
# % keywords: elevation
# % keywords: slopeunits
# %end

# %option G_OPT_R_INPUT
# % key: demmap
# % description: Input passed to r.slopeunits.create and metrics
# % required : yes
# %end

# %option G_OPT_R_INPUT
# % key: plainsmap
# % description: Input passed to r.slopeunits.create
# % required : yes
# %end

# %option
# % key: thresh
# % type: double
# % description: Input passed to r.slopeunits.create
# % required : yes
# %end

# %option
# % key: rf
# % type: integer
# % description: Input passed to r.slopeunits.create
# % required : yes
# %end

# %option
# % key: maxiteration
# % type: integer
# % description: Input passed to r.slopeunits.create
# % required : yes
# %end

# %option
# % key: cleansize
# % type: double
# % answer: 25000
# % description: Input passed to r.slopeunits.create and metrics
# % required: yes
# %end

# %option
# % key: basin
# % type: string
# % description: Input passed to r.slopeunits.metrics
# % required: yes
# %end

import atexit
import math
import os

import grass.script as grass

# TODO make configurable or get rid of all temp files
outdir = "/workdir/test/outdir_test_05_py"

# initialize global vars
rm_rasters = []
rm_vectors = []
count_global = 0


def cleanup():
    """Cleanup fuction"""
    nuldev = open(os.devnull, "w")
    kwargs = {"flags": "f", "quiet": True, "stderr": nuldev}
    for rmrast in rm_rasters:
        if grass.find_file(name=rmrast, element="cell")["file"]:
            grass.run_command("g.remove", type="raster", name=rmrast, **kwargs)
    for rmvect in rm_vectors:
        if grass.find_file(name=rmvect, element="vector")["file"]:
            grass.run_command("g.remove", type="vector", name=rmvect, **kwargs)

    if grass.find_file("MASK")["file"]:
        # grass.run_command('r.mask', flags='r')
        grass.run_command(
            "g.remove", type="raster", name="MASK", flags="f", quiet=True
        )


def calcola_loop(
    x, y, outdir, basin, dem, rf, maxiteration, thresh, cleansize, plainsmap
):
    """_summary_

    Args:
        x (_type_): _description_
        y (_type_): _description_
        outdir (_type_): _description_
    """

    global count_global

    env = grass.gisenv()
    gisdbase = env["GISDBASE"]
    location = env["LOCATION_NAME"]
    master_mapset = env["MAPSET"]

    # we check if this point was calculated already. We only check for V
    with open(os.path.join(outdir, "calcd.dat"), "r") as file:
        found_V = next(
            (
                float(line.split()[2])
                for line in file
                if float(line.split()[0]) == x and float(line.split()[1]) == y
            ),
            None,
        )

    # found_V will be None if no match was found
    if not found_V:
        # Not calculated before. We call run_batch and calculate this point
        count_global += 1
        cval = x
        aval = y
        mapset_prefix = "tmp_su"
        # TODO is this needed?
        ico = count_global
        nome_mapset = f"{mapset_prefix}_{ico}"

        grass.utils.try_rmdir(os.path.join(gisdbase, location, nome_mapset))
        grass.run_command("g.mapset", mapset=nome_mapset, flags="c")
        grass.run_command("g.mapsets", mapset=master_mapset, operation="add")
        grass.run_command(
            "g.region",
            vect=f"{basin}@{master_mapset}",
            align=f"{dem}@{master_mapset}",
        )

        grass.message(
            f"Calculating slopeunits for x={str(x)} and y={str(y)} ..."
        )

        # TODO use GRASS GIS in python with our addons when they are ready
        import subprocess

        process = subprocess.Popen(
            [
                "/workdir/test/r.slopeunits",
                f"demmap={dem}@{master_mapset}",
                "slumap=su_tmp",
                "slumapclean=su_tmp_cl",
                f"rf={rf}",
                f"maxiteration={maxiteration}",
                f"thresh={thresh}",
                f"areamin={aval}",
                f"cvmin={cval}",
                f"cleansize={cleansize}",
                f"plainsmap={plainsmap}@{master_mapset}",
                "-m",
                "--o",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = process.communicate()
        stdout
        # stderr

        # grass.run_command(
        #     "r.slopeunits.create",
        #     demmap=f"{dem}@{master_mapset}",
        #     plainsmap=f"{plainsmap}@{master_mapset}",
        #     slumap="su_tmp",
        #     thresh=thresh,
        #     areamin=aval,
        #     cvmin=cval,
        #     rf=rf,
        #     maxiteration=maxiteration,
        #     overwrite=True,
        # )
        # grass.run_command(
        #     "r.slopeunits.clean",
        #     demmap=f"{dem}@{master_mapset}",
        #     plainsmap=f"{plainsmap}@{master_mapset}",
        #     slumap="su_tmp",
        #     slumapclean="su_tmp_cl",
        #     cleansize=cleansize,
        #     flags=["m"],
        # )

        region = grass.parse_command("g.region", flags="pg")
        resolution = math.floor(
            float(region["ewres"]) * float(region["nsres"])
        )
        outfile = f"{outdir}/objf_{ico}.dat"

        grass.run_command(
            "r.to.vect",
            input="su_tmp_cl",
            output="su_tmp_cl",
            type="area",
            overwrite=True,
        )
        # TODO: clean with addon cleanup function
        grass.run_command(
            "g.remove",
            type="raster",
            name="su_tmp,su_tmp_cl",
            flags="f",
        )

        grass.message(f"Calculating metrics for x={str(x)} and y={str(y)} ...")

        # TODO use GRASS GIS in python with our addon when ready
        import subprocess

        process = subprocess.Popen(
            [
                "python",
                "/workdir/test/calculate_metric_alt2.py",
                f"{basin}",
                f"{dem}",
                "su_tmp_cl",
                f"{cleansize}",
                f"{ico}",
                f"{aval}",
                f"{cval}",
                f"{resolution}",
                f"{outfile}",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = process.communicate()
        stdout
        # stderr

        # grass.run_command(
        #     "r.slopeunits.metrics",
        #     basin=f{basin}",
        #     dem=f"{dem}",
        #     sucl="su_tmp_cl",
        #     thr_clean=f"{cleansize}",
        #     uid=f"{ico}",
        #     areamin=f"{aval}",
        #     cvmin=f"{cval}",
        #     res=f"{resolution}",
        #     outfile=f"{outfile}",
        # )

        grass.message(f"exe no. ${ico}: done")

        grass.run_command("g.mapset", mapset=master_mapset)
        grass.utils.try_rmdir(os.path.join(gisdbase, location, nome_mapset))

        file_path = os.path.join(outdir, f"objf_{count_global}.dat")
        with open(file_path, "r") as file:
            out1 = f"{float(file.read().split()[2]):16.14f}"
        with open(file_path, "r") as file:
            out2 = f"{float(file.read().split()[3]):16.9f}".strip()
        with open(os.path.join(outdir, "tmp.txt"), "a") as file:
            file.write(f"{out1} {out2} {count_global}\n")

        grass.message(f"Writing to calcd.dat: {x} {y} {out1} {out2} ...")
        with open(os.path.join(outdir, "calcd.dat"), "a") as file:
            file.write(f"{x} {y} {out1} {out2}\n")

    else:
        # found_V found - x and y already processed and in file
        grass.message(
            "x and y already calculated! Do nothing but copy the result "
            "in the result vectors"
        )

        out1 = float(f"{found_V:16.14f}")
        with open(os.path.join(outdir, "calcd.dat"), "r") as file:
            found_I = next(
                (
                    float(line.split()[3])
                    for line in file
                    if float(line.split()[0]) == x
                    and float(line.split()[1]) == y
                ),
                None,
            )
        out2 = float(f"{found_I:16.14f}")

    grass.message(f"Writing to current.txt: {ico} {x} {y} {out1} {out2} ...")
    with open(os.path.join(outdir, "current.txt"), "a") as file:
        file.write(f"{ico} {x} {y} {out1} {out2}\n")


def calcola_current(
    x_lims,
    y_lims,
    outdir,
    basin,
    dem,
    rf,
    maxiteration,
    thresh,
    cleansize,
    plainsmap,
):
    """_summary_

    Args:
        x_lims (_type_): _description_
        y_lims (_type_): _description_
        outdir (_type_): _description_
    """

    ico = 0

    # TODO check if they need to be duplicated
    x_lims_local = x_lims
    y_lims_local = y_lims

    for x in x_lims_local:
        for y in y_lims_local:
            # x = 0.05000000000000
            # y = 50000.0000000000
            calcola_loop(
                x,
                y,
                outdir,
                basin,
                dem,
                rf,
                maxiteration,
                thresh,
                cleansize,
                plainsmap,
            )
            ico += 1

    # now the fifth, central point

    # # TODO: check if precision is important !?
    # import decimal
    # # Set precision for decimal calculations
    # decimal.getcontext().prec = 14
    # # Assuming x_lims_local and y_lims_local are lists with two elements each
    # x_tmp = (
    # decimal.Decimal(x_lims_local[0]) + decimal.Decimal(x_lims_local[1])) / 2
    # y_tmp = (
    # decimal.Decimal(y_lims_local[0]) + decimal.Decimal(y_lims_local[1])) / 2
    # x_half = f"{x_tmp:.14f}"
    # y_half = f"{y_tmp:.9f}"
    # # 0.15000000000000 125000 0.15000000000000 125000.000000000

    x_tmp = (x_lims_local[0] + x_lims_local[1]) / 2
    y_tmp = (y_lims_local[0] + y_lims_local[1]) / 2
    x_half = f"{x_tmp:16.14f}"
    y_half = f"{y_tmp:16.9f}"

    calcola_loop(
        x_half,
        y_half,
        outdir,
        basin,
        dem,
        rf,
        maxiteration,
        thresh,
        cleansize,
        plainsmap,
    )


def main():
    """Main function of r.slopeunits.metrics"""
    global rm_rasters, rm_vectors
    global count_global

    dem = options["demmap"]
    plainsmap = options["plainsmap"]
    thresh = float(options["thresh"])
    rf = int(options["rf"])
    maxiteration = int(options["maxiteration"])
    cleansize = options["cleansize"]
    basin = options["basin"]

    # Clean start
    grass.utils.try_rmdir(outdir)
    os.mkdir(outdir)
    with open(os.path.join(outdir, "calcd.dat"), "w") as file:
        file.write("")

    # Start search: initial limits (x is cvar; y is amin)
    x_lims = [0.05000000000000, 0.25000000000000]
    y_lims = [50000.0000000000, 200000.000000000]

    # TODO: check if setting these is needed
    x_lims_old = x_lims
    y_lims_old = y_lims
    x_lims_cur = x_lims_old
    y_lims_cur = y_lims_old

    epsilonx = 0.01
    epsilony = 50000
    istop = 0
    count = 0

    while istop == 0:
        # this function calculates slope units & corresponding F(a,c)
        # one the four current points
        calcola_current(
            x_lims,
            y_lims,
            outdir,
            basin,
            dem,
            rf,
            maxiteration,
            thresh,
            cleansize,
            plainsmap,
        )
        count += 1

        # TODO continue translating test05.sh from line 184
        # now we need to find V_min,V_max,I_min,I_max corresponding ...
        # ...
        # Requirements to continue:
        # - calcd.dat
        # - current.txt


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
