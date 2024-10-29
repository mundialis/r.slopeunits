#!/usr/bin/env python3
#
############################################################################
#
# MODULE:       r.slopeunits.optimize for GRASS 8
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
# % description: Optimize inputs for slope units
# % keywords: raster
# % keywords: elevation
# % keywords: slopeunits
# %end

# %option G_OPT_R_INPUT
# % key: demmap
# % description: Input digital elevation model
# % required : yes
# %end

# %option G_OPT_R_INPUT
# % key: plainsmap
# % description: Input raster map of alluvial plains
# % required : yes
# %end

# %option G_OPT_V_INPUT
# % key: basin
# % description: Input basin passed to r.slopeunits.metrics
# % required: yes
# %end

# %option G_OPT_R_OUTPUT
# % key: slumap
# % description: Slope Units layer (the main output of r.slopeunits.create)
# % required : yes
# % answer: su_tmp
# %end

# %option G_OPT_R_OUTPUT
# % key: slumapclean
# % description: Slope Units layer, cleaned (the main output of r.slopeunits.clean)
# % required : no
# % answer: su_tmp_cl
# %end

# %option
# % key: thresh
# % type: double
# % description: Initial threshold (m^2)
# % required : yes
# %end

# %option
# % key: rf
# % type: integer
# % description: Factor used to iterativelly reduce initial threshold
# % required : yes
# %end

# %option
# % key: maxiteration
# % type: integer
# % description: maximum number of iteration to do before stopped
# % required : yes
# %end

# %option
# % key: cleansize
# % type: double
# % answer: 25000
# % description: Slope Units size to be removed
# % required: yes
# %end

# %option
# % key: cvmin
# % type: double
# % answer: 0.05,0.25
# % description: Start search with these initial minimum and maximum values of the circular variance (0.0-1.0) below which the slope unit is not further segmented
# % multiple : no
# % key_desc : min,max
# % required : yes
# %end

# %option
# % key: areamin
# % type: double
# % answer: 50000.0,200000.0
# % description: Start search with these initial minimum and maximum values of the area (m^2) below which the slope unit is not further segmented
# % multiple : no
# % key_desc : min,max
# % required : yes
# %end

# %option
# % key: epsilonx
# % type: double
# % answer: 0.01
# % description: Stop loop when difference of cvmin limits is lesser than this value
# % required: yes
# %end

# %option
# % key: epsilony
# % type: double
# % answer: 50000
# % description: Stop loop when difference of areamin limits is lesser than this value
# % required: yes
# %end

# %option G_OPT_M_DIR
# % key: outdir
# % description: Output directory for intermediate results for all cvmin and areamin values and final results. Default folder "outdir" in current working directory.
# % required: yes
# % answer: outdir
# %end

# pylint: disable=C0302 (too-many-lines)

import atexit
import math
import os

import grass.script as grass

# initialize global vars
rm_rasters = []
rm_vectors = []
COUNT_GLOBAL = 0
env = grass.gisenv()
gisdbase = env["GISDBASE"]
location = env["LOCATION_NAME"]
master_mapset = env["MAPSET"]


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


def run_batch(
    cvmin,
    areamin,
    basin,
    dem,
    slumap,
    slumapclean,
    redf,
    maxiteration,
    thresh,
    cleansize,
    plainsmap,
):
    """Calls r.slopeunits.create, r.slopeunits.clean and r.slopeunits.metrics"""
    global COUNT_GLOBAL
    global gisdbase
    global location
    global master_mapset

    COUNT_GLOBAL += 1
    ico = COUNT_GLOBAL
    mapset_prefix = "tmp_su"
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
        f"Calculating slopeunits for cvmin={str(cvmin)} and "
        f"areamin={str(areamin)} ..."
    )
    grass.run_command(
        "r.slopeunits.create",
        demmap=f"{dem}@{master_mapset}",
        plainsmap=f"{plainsmap}@{master_mapset}",
        slumap=slumap,
        thresh=thresh,
        areamin=areamin,
        cvmin=cvmin,
        rf=redf,
        maxiteration=maxiteration,
        overwrite=True,
    )
    grass.run_command(
        "r.slopeunits.clean",
        demmap=f"{dem}@{master_mapset}",
        plainsmap=f"{plainsmap}@{master_mapset}",
        slumap=slumap,
        slumapclean=slumapclean,
        cleansize=cleansize,
        flags=["m"],
        overwrite=True,
    )

    region = grass.parse_command("g.region", flags="pg")
    resolution = math.floor(float(region["ewres"]) * float(region["nsres"]))

    grass.run_command(
        "r.to.vect",
        input=slumapclean,
        output=slumapclean,
        type="area",
        overwrite=True,
        quiet=True,
    )
    grass.run_command(
        "g.remove",
        type="raster",
        name=f"{slumap},{slumapclean}",
        flags="f",
        quiet=True,
    )

    grass.message(
        f"Calculating metrics for cvmin={str(cvmin)} and "
        "areamin={str(areamin)} ..."
    )
    metrics = grass.run_command(
        "r.slopeunits.metrics",
        basin=basin,
        dem=dem,
        slumapclean=slumapclean,
        cleansize=cleansize,
        areamin=areamin,
        cvmin=cvmin,
        resolution=resolution,
    )

    grass.message(f"exe no. {ico}: done")

    grass.run_command("g.mapset", mapset=master_mapset)
    grass.utils.try_rmdir(os.path.join(gisdbase, location, nome_mapset))

    return metrics


def calcola_loop(
    cvmin,
    areamin,
    basin,
    dem,
    slumap,
    slumapclean,
    redf,
    maxiteration,
    thresh,
    cleansize,
    plainsmap,
    ico,
    calcd_file,
    current_file,
):
    """Wrapper for calculating slopeunits and metrics to parse the data"""

    global COUNT_GLOBAL

    # we check if this point was calculated already. We only check for V
    with open(calcd_file, "r") as file:
        found_v = next(
            (
                float(line.split()[2])
                for line in file
                if float(line.split()[0]) == cvmin
                and float(line.split()[1]) == areamin
            ),
            None,
        )

    # found_v will be None if no match was found
    if not found_v:
        # Not calculated before. We call run_batch and calculate this point
        metrics = run_batch(
            cvmin,
            areamin,
            basin,
            dem,
            slumap,
            slumapclean,
            redf,
            maxiteration,
            thresh,
            cleansize,
            plainsmap,
        )
        out1 = f"{float(metrics['v_fin']):16.14f}"
        out2 = f"{float(metrics['i_fin']):16.9f}".strip()
        grass.message(
            f"Writing to calcd.dat: {cvmin} {areamin} {out1} {out2} ..."
        )
        with open(calcd_file, "a") as file:
            file.write(f"{cvmin} {areamin} {out1} {out2}\n")

    else:
        # found_v found - x and y already processed and in file
        grass.message(
            "x and y already calculated! Do nothing but copy the result "
            "in the result vectors"
        )

        out1 = float(f"{found_v:16.14f}")
        with open(calcd_file, "r") as file:
            found_i = next(
                (
                    float(line.split()[3])
                    for line in file
                    if float(line.split()[0]) == cvmin
                    and float(line.split()[1]) == areamin
                ),
                None,
            )
        out2 = float(f"{found_i:16.14f}")

    grass.message(
        f"Writing to current.txt: {ico} {cvmin} {areamin} {out1} {out2} ..."
    )
    with open(current_file, "a") as file:
        file.write(f"{ico} {cvmin} {areamin} {out1} {out2}\n")


def calcola_current(
    x_lims,
    y_lims,
    basin,
    dem,
    slumap,
    slumapclean,
    redf,
    maxiteration,
    thresh,
    cleansize,
    plainsmap,
    calcd_file,
    current_file,
):
    """Calculates slope units & corresponding F(a,c)
    on the four current points and the fifth, central point
    """

    ico = 0

    x_lims_local = x_lims.copy()
    y_lims_local = y_lims.copy()

    for cvmin in x_lims_local:
        for areamin in y_lims_local:
            calcola_loop(
                cvmin,
                areamin,
                basin,
                dem,
                slumap,
                slumapclean,
                redf,
                maxiteration,
                thresh,
                cleansize,
                plainsmap,
                ico,
                calcd_file,
                current_file,
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

    x_half = (x_lims_local[0] + x_lims_local[1]) / 2
    y_half = (y_lims_local[0] + y_lims_local[1]) / 2
    # x_half = f"{x_tmp:16.14f}"
    # y_half = f"{y_tmp:16.9f}"
    calcola_loop(
        x_half,
        y_half,
        basin,
        dem,
        slumap,
        slumapclean,
        redf,
        maxiteration,
        thresh,
        cleansize,
        plainsmap,
        ico,
        calcd_file,
        current_file,
    )


def calculate_new_limits_caso_1(im1, im2, x_lims_cur, y_lims_cur):
    """_summary_

    Args:
        im1 (_type_): _description_
        im2 (_type_): _description_
        x_lims_cur (_type_): _description_
        y_lims_cur (_type_): _description_

    Returns:
        _type_: _description_
    """
    caso = ""
    x_lims_new = [0, 0]
    y_lims_new = [0, 0]

    # my first approach was to use a rectangle on the four points
    #     x_lims_new[0]=$x_min
    #     x_lims_new[1]=$x_max
    #     y_lims_new[0]=$y_min
    #     y_lims_new[1]=$y_max
    # now I use a slighly larger rectangle
    if im1 == 2 or im2 == 2:
        caso = "1a)"
        tmp1 = (2.0 * x_lims_cur[0] + x_lims_cur[1]) / 3.0
        x_lims_new[0] = float(f"{tmp1:.14f}")
        x_lims_new[1] = x_lims_cur[1]
        y_lims_new[0] = y_lims_cur[0]
        tmp1 = (y_lims_cur[0] + 2.0 * y_lims_cur[1]) / 3.0
        y_lims_new[1] = float(f"{tmp1:.9f}")

    elif im1 == 3 or im2 == 3:
        caso = "1b)"
        tmp1 = (2.0 * x_lims_cur[0] + x_lims_cur[1]) / 3.0
        x_lims_new[0] = float(f"{tmp1:.14f}")
        x_lims_new[1] = x_lims_cur[1]
        tmp1 = (2.0 * y_lims_cur[0] + y_lims_cur[1]) / 3.0
        y_lims_new[0] = float(f"{tmp1:.9f}")
        y_lims_new[1] = y_lims_cur[1]

    elif im1 == 1 or im2 == 1:
        caso = "1c)"
        x_lims_new[0] = x_lims_cur[0]
        tmp1 = (x_lims_cur[0] + 2.0 * x_lims_cur[1]) / 3.0
        x_lims_new[1] = float(f"{tmp1:.14f}")
        tmp1 = (2.0 * y_lims_cur[0] + y_lims_cur[1]) / 3.0
        y_lims_new[0] = float(f"{tmp1:.9f}")
        y_lims_new[1] = y_lims_cur[1]

    elif im1 == 0 or im2 == 0:
        caso = "1d)"
        x_lims_new[0] = x_lims_cur[0]
        tmp1 = (x_lims_cur[0] + 2.0 * x_lims_cur[1]) / 3.0
        x_lims_new[1] = float(f"{tmp1:.14f}")
        y_lims_new[0] = y_lims_cur[0]
        tmp1 = (y_lims_cur[0] + 2.0 * y_lims_cur[1]) / 3.0
        y_lims_new[1] = float(f"{tmp1:.9f}")

    return caso, x_lims_new, y_lims_new


def calculate_new_limits_caso_2(yp1, x_lims_cur, y_lims_cur):
    """_summary_

    Args:
        yp1 (_type_): _description_
        x_lims_cur (_type_): _description_
        y_lims_cur (_type_): _description_

    Returns:
        _type_: _description_
    """
    caso = ""
    x_lims_new = [0, 0]
    y_lims_new = [0, 0]

    deltay1 = yp1 - y_lims_cur[0]
    deltay2 = deltay1 * deltay1
    deltay1 = math.sqrt(deltay2)

    if deltay1 > 0:
        # y^1=y^2=y_1=y_3: caso 2b)
        caso = "2b)"
        tmp1 = (y_lims_cur[0] + 2.0 * y_lims_cur[1]) / 3.0
        y_lims_new[0] = float(f"{tmp1:.9f}")
        y_lims_new[1] = y_lims_cur[1]
    else:
        # y^1=y^2=y_0=y_2: caso 2a)
        caso = "2a)"
        y_lims_new[0] = y_lims_cur[0]
        tmp1 = (2.0 * y_lims_cur[0] + y_lims_cur[1]) / 3.0
        y_lims_new[1] = float(f"{tmp1:.9f}")

    x_lims_new[0] = x_lims_cur[0]
    x_lims_new[1] = x_lims_cur[1]

    return caso, x_lims_new, y_lims_new


def calculate_new_limits_caso_3(xp1, x_lims_cur, y_lims_cur):
    """_summary_

    Args:
        xp1 (_type_): _description_
        x_lims_cur (_type_): _description_
        y_lims_cur (_type_): _description_

    Returns:
        _type_: _description_
    """
    caso = ""
    x_lims_new = [0, 0]
    y_lims_new = [0, 0]

    deltax1 = xp1 - x_lims_cur[0]
    deltax2 = deltax1 * deltax1
    deltax1 = math.sqrt(deltax2)

    if deltax1 > 0:
        # x^1=x^2=x_2=x_3: caso 3b)
        caso = "3b)"
        tmp1 = (2.0 * x_lims_cur[0] + x_lims_cur[1]) / 3.0
        x_lims_new[0] = float(f"{tmp1:.14f}")
        x_lims_new[1] = x_lims_cur[1]
    else:
        # x^1=x^2=x_0=x_1: caso 3a)
        caso = "3a)"
        x_lims_new[0] = x_lims_cur[0]
        tmp1 = (x_lims_cur[0] + 2.0 * x_lims_cur[1]) / 3.0
        x_lims_new[1] = float(f"{tmp1:.14f}")

    y_lims_new[0] = y_lims_cur[0]
    y_lims_new[1] = y_lims_cur[1]

    return caso, x_lims_new, y_lims_new


def calculate_new_limits_caso_4(x_lims_cur, y_lims_cur):
    """_summary_

    Args:
        x_lims_cur (_type_): _description_
        y_lims_cur (_type_): _description_

    Returns:
        _type_: _description_
    """
    caso = "4)"
    # Copy the current limits to new limits
    x_lims_new = x_lims_cur.copy()
    y_lims_new = y_lims_cur.copy()

    return caso, x_lims_new, y_lims_new


def main():
    """Main function of r.slopeunits.optimize"""
    global rm_rasters, rm_vectors
    global COUNT_GLOBAL

    dem = options["demmap"]
    plainsmap = options["plainsmap"]
    slumap = options["slumap"]
    slumapclean = options["slumapclean"]
    thresh = float(options["thresh"])
    redf = int(options["rf"])
    maxiteration = int(options["maxiteration"])
    cleansize = options["cleansize"]
    basin = options["basin"]
    x_lims = [
        float(options["cvmin"].split(",")[0]),
        float(options["cvmin"].split(",")[1]),
    ]
    y_lims = [
        float(options["areamin"].split(",")[0]),
        float(options["areamin"].split(",")[1]),
    ]
    epsilonx = options["epsilonx"]
    epsilony = options["epsilony"]
    outdir = os.path.abspath(options["outdir"])

    calcd_file = os.path.join(outdir, "calcd.dat")
    current_file = os.path.join(outdir, "current.txt")
    optimum_file = os.path.join(outdir, "opt.txt")

    # Clean start
    if os.path.exists(outdir):
        if grass.overwrite():
            grass.utils.try_remove(calcd_file)
            grass.utils.try_remove(current_file)
            grass.utils.try_remove(optimum_file)
        else:
            if (
                os.path.exists(calcd_file)
                or os.path.exists(current_file)
                or os.path.exists(optimum_file)
            ):
                grass.fatal(
                    f"ERROR: One of {calcd_file}, {current_file} or "
                    f"{optimum_file} exists. To overwrite, use the --overwrite "
                    "flag"
                )
    else:
        os.mkdir(outdir)

    # File will be read before written, so empty file need to exist
    with open(calcd_file, "w") as file:
        file.write("")

    x_lims_cur = x_lims.copy()
    y_lims_cur = y_lims.copy()

    istop = 0
    count = 0

    while istop == 0:
        # clean current_file before method call
        with open(current_file, "w"):
            pass

        calcola_current(
            x_lims_cur,
            y_lims_cur,
            basin,
            dem,
            slumap,
            slumapclean,
            redf,
            maxiteration,
            thresh,
            cleansize,
            plainsmap,
            calcd_file,
            current_file,
        )
        count += 1

        # now we need to find v_min,v_max,i_min,i_max corresponding ...
        v_min, v_max, i_min, i_max = (
            float("inf"),
            float("-inf"),
            float("inf"),
            float("-inf"),
        )
        with open(calcd_file, "r") as file:
            for line in file:
                fields = line.strip().split()
                if len(fields) >= 4:
                    v_cur = float(fields[2])
                    i_cur = float(fields[3])
                    v_min = min(v_min, v_cur)
                    v_max = max(v_max, v_cur)
                    i_min = min(i_min, i_cur)
                    i_max = max(i_max, i_cur)

        # ... and use them to calculate the metric for the 5 punti of the last
        # calculate rectangle
        results = []
        with open(current_file, "r") as file:
            for line in file:
                fields = line.strip().split()
                if len(fields) >= 5:
                    v_cur = float(fields[3])
                    i_cur = float(fields[4])
                    calculated_value = (v_max - v_cur) / (v_max - v_min) + (
                        (i_max - i_cur) / (i_max - i_min)
                    )
                    results.append(fields + [str(calculated_value)])
        # Sort the results based on the calculated value (6th column)
        sorted_results = sorted(results, key=lambda x: str(x[5]))
        with open(current_file, "w") as file:
            for result in sorted_results:
                file.write(" ".join(result) + "\n")

        im1 = float(sorted_results[-1][0])
        xp1 = float(sorted_results[-1][1])
        yp1 = float(sorted_results[-1][2])
        im2 = float(sorted_results[-2][0])
        xp2 = float(sorted_results[-2][1])
        yp2 = float(sorted_results[-2][2])

        deltax1 = float(xp1) - float(xp2)
        deltax2 = deltax1 * deltax1
        deltax = math.sqrt(deltax2)
        deltay1 = float(yp1) - float(yp2)
        deltay2 = deltay1 * deltay1
        deltay = math.sqrt(deltay2)
        # Format the result to 12 decimal places
        # formatted_deltax = "{:.12f}".format(deltax)

        # see comment to first approach in calculate_new_limits_caso_1
        # x_min = xp1 if xp1 < xp2 else xp2
        # x_max = xp1 if xp1 > xp2 else xp2
        # y_min = yp1 if yp1 < yp2 else yp2
        # y_max = yp1 if yp1 > yp2 else yp2
        print(im1, im2, x_lims_cur, y_lims_cur, xp1, xp2, yp1, yp2)

        # we use deltax to distinguish if x1=x2 or x1!=x2; same for y
        x_lims_new = [0, 0]
        y_lims_new = [0, 0]
        caso = ""
        if deltax > 0:
            if deltay > 0:
                # both different
                print(count, "(XdYd)")
                caso, x_lims_new, y_lims_new = calculate_new_limits_caso_1(
                    im1, im2, x_lims_cur, y_lims_cur
                )
            else:
                # x different, y same
                print(count, "(XdYe)")
                caso, x_lims_new, y_lims_new = calculate_new_limits_caso_2(
                    yp1, x_lims_cur, y_lims_cur
                )
        else:
            if deltay > 0:
                # x same, y different
                print(count, "(XeYd)")
                caso, x_lims_new, y_lims_new = calculate_new_limits_caso_3(
                    xp1, x_lims_cur, y_lims_cur
                )
            else:
                # both same
                print(count, "(XeYe)")
                caso, x_lims_new, y_lims_new = calculate_new_limits_caso_4(
                    x_lims_cur, y_lims_cur
                )
                istop = 1

        print(
            f"{caso} - {str(x_lims_new[0])} {str(x_lims_new[1])} "
            f"{str(y_lims_new[0])} {str(y_lims_new[1])}"
        )
        # with open(os.path.join(outdir, "log.txt"), "a") as file:
        #     file.write(
        #         f"{caso} - {str(x_lims_new[0])} {str(x_lims_new[1])} "
        #         f"{str(y_lims_new[0])} {str(y_lims_new[1])}\n"
        #     )

        if istop >= 0:
            x_lims_cur = x_lims_new.copy()
            y_lims_cur = y_lims_new.copy()

            deltax1 = x_lims_new[0] - x_lims_new[1]
            deltax2 = deltax1 * deltax1
            deltax = math.sqrt(deltax2)

            deltay1 = y_lims_new[0] - y_lims_new[1]
            deltay2 = deltay1 * deltay1
            deltay = math.sqrt(deltay2)
            if deltax < epsilonx:
                if deltay < epsilony:
                    istop = 1

    # end while
    x_opt = (x_lims_cur[0] + x_lims_cur[1]) / 2.0
    y_opt = (y_lims_cur[0] + y_lims_cur[1]) / 2.0
    print(f"x_opt: {x_opt}")
    print(f"y_opt: {y_opt}")

    with open(optimum_file, "w") as file:
        file.write(f"x_opt: {x_opt}\n")
        file.write(f"y_opt: {y_opt}\n")

    with open(current_file, "r") as file:
        print(file.read())

    # Print function evaluations
    print(f" valutazioni funzione: {COUNT_GLOBAL}")
    # End of search


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
