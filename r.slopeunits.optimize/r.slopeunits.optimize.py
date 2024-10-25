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


import atexit
import math
import os

import grass.script as grass

# TODO make configurable or get rid of all temp files
outdir = "/workdir/test/outdir_test05_py"
# TODO make configurable?
#         "slumap=su_tmp",
#         "slumapclean=su_tmp_cl",

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
    x, y, outdir, basin, dem, redf, maxiteration, thresh, cleansize, plainsmap
):
    """Calls r.slopeunits.create, r.slopeunits.clean and r.slopeunits.metrics

    Args:
        x (_type_): _description_
        y (_type_): _description_
        outdir (_type_): _description_
        basin (_type_): _description_
        dem (_type_): _description_
        redf (_type_): _description_
        maxiteration (_type_): _description_
        thresh (_type_): _description_
        cleansize (_type_): _description_
        plainsmap (_type_): _description_
    """
    global COUNT_GLOBAL
    global gisdbase
    global location
    global master_mapset

    COUNT_GLOBAL += 1
    # cval = x
    # aval = y
    mapset_prefix = "tmp_su"
    ico = COUNT_GLOBAL
    nome_mapset = f"{mapset_prefix}_{ico}"

    grass.utils.try_rmdir(os.path.join(gisdbase, location, nome_mapset))
    grass.run_command("g.mapset", mapset=nome_mapset, flags="c")
    grass.run_command("g.mapsets", mapset=master_mapset, operation="add")
    grass.run_command(
        "g.region",
        vect=f"{basin}@{master_mapset}",
        align=f"{dem}@{master_mapset}",
    )

    grass.message(f"Calculating slopeunits for x={str(x)} and y={str(y)} ...")

    grass.run_command(
        "r.slopeunits.create",
        demmap=f"{dem}@{master_mapset}",
        plainsmap=f"{plainsmap}@{master_mapset}",
        slumap="su_tmp",
        thresh=thresh,
        areamin=y,
        cvmin=x,
        rf=redf,
        maxiteration=maxiteration,
        overwrite=True,
    )
    grass.run_command(
        "r.slopeunits.clean",
        demmap=f"{dem}@{master_mapset}",
        plainsmap=f"{plainsmap}@{master_mapset}",
        slumap="su_tmp",
        slumapclean="su_tmp_cl",
        cleansize=cleansize,
        flags=["m"],
        overwrite=True,
    )

    region = grass.parse_command("g.region", flags="pg")
    resolution = math.floor(float(region["ewres"]) * float(region["nsres"]))
    outfile = f"{outdir}/objf_{ico}.dat"

    grass.run_command(
        "r.to.vect",
        input="su_tmp_cl",
        output="su_tmp_cl",
        type="area",
        overwrite=True,
        quiet=True,
    )
    grass.run_command(
        "g.remove",
        type="raster",
        name="su_tmp,su_tmp_cl",
        flags="f",
        quiet=True,
    )

    grass.message(f"Calculating metrics for x={str(x)} and y={str(y)} ...")
    grass.run_command(
        "r.slopeunits.metrics",
        basin=basin,
        dem=dem,
        slumapclean="su_tmp_cl",
        cleansize=cleansize,
        areamin=y,
        cvmin=x,
        resolution=resolution,
        outfile=outfile,
    )

    grass.message(f"exe no. {ico}: done")

    grass.run_command("g.mapset", mapset=master_mapset)
    grass.utils.try_rmdir(os.path.join(gisdbase, location, nome_mapset))


def calcola_loop(
    x,
    y,
    outdir,
    basin,
    dem,
    redf,
    maxiteration,
    thresh,
    cleansize,
    plainsmap,
    ico,
):
    """_summary_

    Args:
        x (_type_): _description_
        y (_type_): _description_
        outdir (_type_): _description_
        basin (_type_): _description_
        dem (_type_): _description_
        redf (_type_): _description_
        maxiteration (_type_): _description_
        thresh (_type_): _description_
        cleansize (_type_): _description_
        plainsmap (_type_): _description_
        ico (_type_): _description_
    """

    global COUNT_GLOBAL

    # we check if this point was calculated already. We only check for V
    with open(os.path.join(outdir, "calcd.dat"), "r") as file:
        found_v = next(
            (
                float(line.split()[2])
                for line in file
                if float(line.split()[0]) == x and float(line.split()[1]) == y
            ),
            None,
        )

    # found_v will be None if no match was found
    if not found_v:
        # Not calculated before. We call run_batch and calculate this point
        run_batch(
            x,
            y,
            outdir,
            basin,
            dem,
            redf,
            maxiteration,
            thresh,
            cleansize,
            plainsmap,
        )
        file_path = os.path.join(outdir, f"objf_{COUNT_GLOBAL}.dat")
        with open(file_path, "r") as file:
            out1 = f"{float(file.read().split()[2]):16.14f}"
        with open(file_path, "r") as file:
            out2 = f"{float(file.read().split()[3]):16.9f}".strip()

        grass.message(f"Writing to calcd.dat: {x} {y} {out1} {out2} ...")
        with open(os.path.join(outdir, "calcd.dat"), "a") as file:
            file.write(f"{x} {y} {out1} {out2}\n")

    else:
        # found_v found - x and y already processed and in file
        grass.message(
            "x and y already calculated! Do nothing but copy the result "
            "in the result vectors"
        )

        out1 = float(f"{found_v:16.14f}")
        with open(os.path.join(outdir, "calcd.dat"), "r") as file:
            found_i = next(
                (
                    float(line.split()[3])
                    for line in file
                    if float(line.split()[0]) == x
                    and float(line.split()[1]) == y
                ),
                None,
            )
        out2 = float(f"{found_i:16.14f}")

    grass.message(f"Writing to current.txt: {ico} {x} {y} {out1} {out2} ...")
    with open(os.path.join(outdir, "current.txt"), "a") as file:
        file.write(f"{ico} {x} {y} {out1} {out2}\n")


def calcola_current(
    x_lims,
    y_lims,
    outdir,
    basin,
    dem,
    redf,
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
        basin (_type_): _description_
        dem (_type_): _description_
        redf (_type_): _description_
        maxiteration (_type_): _description_
        thresh (_type_): _description_
        cleansize (_type_): _description_
        plainsmap (_type_): _description_
    """

    ico = 0

    # TODO check if they need to be duplicated
    x_lims_local = x_lims.copy()
    y_lims_local = y_lims.copy()

    for x in x_lims_local:
        for y in y_lims_local:
            calcola_loop(
                x,
                y,
                outdir,
                basin,
                dem,
                redf,
                maxiteration,
                thresh,
                cleansize,
                plainsmap,
                ico,
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

    x_half = (x_lims_local[0] + x_lims_local[1]) / 2
    y_half = (y_lims_local[0] + y_lims_local[1]) / 2
    # x_half = f"{x_tmp:16.14f}"
    # y_half = f"{y_tmp:16.9f}"
    calcola_loop(
        x_half,
        y_half,
        outdir,
        basin,
        dem,
        redf,
        maxiteration,
        thresh,
        cleansize,
        plainsmap,
        ico,
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
    """Main function of r.slopeunits.metrics"""
    global rm_rasters, rm_vectors
    global COUNT_GLOBAL

    dem = options["demmap"]
    plainsmap = options["plainsmap"]
    thresh = float(options["thresh"])
    redf = int(options["rf"])
    maxiteration = int(options["maxiteration"])
    cleansize = options["cleansize"]
    basin = options["basin"]

    # Clean start
    grass.utils.try_rmdir(outdir)
    os.mkdir(outdir)
    # File will be read before written, so empty file need to exist
    with open(os.path.join(outdir, "calcd.dat"), "w") as file:
        file.write("")

    # Start search: initial limits (x is cvar; y is amin)
    x_lims = [0.05000000000000, 0.25000000000000]
    y_lims = [50000.0000000000, 200000.000000000]

    # TODO: check if setting these is needed
    x_lims_old = x_lims.copy()
    y_lims_old = y_lims.copy()
    x_lims_cur = x_lims_old.copy()
    y_lims_cur = y_lims_old.copy()

    epsilonx = 0.01
    epsilony = 50000
    istop = 0
    count = 0

    while istop == 0:
        # clean current.txt before method call
        with open(os.path.join(outdir, "current.txt"), "w"):
            pass

        # this function calculates slope units & corresponding F(a,c)
        # one the four current points
        calcola_current(
            x_lims_cur,
            y_lims_cur,
            outdir,
            basin,
            dem,
            redf,
            maxiteration,
            thresh,
            cleansize,
            plainsmap,
        )
        count += 1

        # now we need to find v_min,v_max,i_min,i_max corresponding ...
        v_min, v_max, i_min, i_max = (
            float("inf"),
            float("-inf"),
            float("inf"),
            float("-inf"),
        )
        with open(os.path.join(outdir, "calcd.dat"), "r") as file:
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
        with open(os.path.join(outdir, "current.txt"), "r") as file:
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

    with open(f"{outdir}/opt.txt", "w") as file:
        file.write(f"x_opt: {x_opt}\n")
        file.write(f"y_opt: {y_opt}\n")

    # TODO remove current.txt? Or write sorted_results to it?
    # with open(f"{outdir}/current.txt", "r") as f:
    #     print(f.read())
    for result in sorted_results:
        print(" ".join(result))

    def print_table(table):
        longest_cols = [
            (max([len(str(row[i])) for row in table]) + 3)
            for i in range(len(table[0]))
        ]
        row_format = "".join(
            ["{:>" + str(longest_col) + "}" for longest_col in longest_cols]
        )
        for row in table:
            print(row_format.format(*row))

    print_table(sorted_results)

    # ico = int(ico) + 1

    # Print function evaluations
    print(f" valutazioni funzione: {COUNT_GLOBAL}")
    # End of search


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
