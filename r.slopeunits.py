#!/usr/bin/env python3
#
############################################################################
#
# MODULE:       r.slopeunits for GRASS 7
# AUTHOR(S):    Ivan Marchesini, Massimiliano Alvioli
# PURPOSE:      To create a raster layer of slope units
# COPYRIGHT:    (C) 2004-2012 by the GRASS Development Team
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
#############################################################################
#
# %module
# % description: Create a raster layer of slope units
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
# % required : no
# %end

# %option G_OPT_R_OUTPUT
# % key: slumap
# % description: Output Slope Units layer (the main output)
# % required : yes
# %end

# %option G_OPT_R_OUTPUT
# % key: slumapclean
# % description: Output Slope Units layer, cleaned (the main output)
# % required : no
# %end

# %option G_OPT_R_OUTPUT
# % key: circvarmap
# % description: Output Circular Variance layer
# % required : no
# %end

# %option G_OPT_R_OUTPUT
# % key: areamap
# % description: Output Area layer; values in square meters
# % required : no
# %end

# %option
# % key: thresh
# % type: double
# % description: Initial threshold (m^2).
# % required : yes
# %end

# %option
# % key: areamin
# % type: double
# % description: Minimum area (m^2) below whitch the slope unit is not further segmented
# % required : yes
# %end

# %option
# % key: areamax
# % type: double
# % description: Maximum area (m^2) above which the slope unit is segmented irrespective of aspect
# % required : no
# %end

# %option
# % key: cvmin
# % type: double
# % description: Minimum value of the circular variance (0.0-1.0) below whitch the slope unit is not further segmented
# % required : yes
# %end

# %option
# % key: rf
# % type: integer
# % description: Factor used to iterativelly reduce initial threshold: newthresh=thresh-thresh/reductionfactor
# % required : yes
# %end

# %option
# % key: maxiteration
# % type: integer
# % description: maximum number of iteration to do before the procedure is in any case stopped
# % required : yes
# %end

# %option
# % key: cleansize
# % type: double
# % description: Slope Units size to be removed
# % required : no
# %end

# %flag
# % key: m
# % description: Perform quick cleaning of small-sized areas and stripes
# % guisection: flags
# %end

# %flag
# % key: n
# % description: Perform detailed cleaning of small-sized areas (slow)
# % guisection: flags
# %end

# import sys
import grass.script as grass
import atexit
import os

# initialize global vars
rm_rasters = []
rm_vectors = []


def cleanup():
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


def slope_units(
    dem=None,
    plains=None,
    slumap=None,
    circvarmap=None,
    areamap=None,
    th=-1.0,
    amin=-1.0,
    cvarmin=-1.0,
    red=-1,
    maxiter=0,
):
    global rm_rasters, rm_vectors
    global thc

    # estimating values of parameters in cells units
    region = grass.region()
    nsres = region["nsres"]
    ewres = region["ewres"]
    thc = int(th / (nsres * ewres))
    aminc = int(amin / (nsres * ewres))
    if options["areamax"]:
        # in square meters
        amax = float(options["areamax"])
        # in cells units
        amaxc = int(amax / (nsres * ewres))

    # a MASK must not exist already
    if grass.find_file(name=rmrast, element="cell")["file"]:
        grass.fatal("Please remove the raster MASK first, before running this module")

    # setting the mask on the DTM
    exp = "$out = if(isnull($mask), null(), 1)"
    grass.mapcalc(exp, out="MASK", mask=dem, quiet=True)
    rm_rasters.append("MASK")

    # generating the aspect layer
    grass.run_command(
        "r.slope.aspect",
        elevation=dem,
        aspect="aspectslutmp",
        quiet=True,
    )
    rm_rasters.append("aspectslutmp")

    # printing something
    grass.message(("Initial threshold (cells) is : %s") % thc)
    grass.message(("Initial minimum area (cells) is : %s") % aminc)

    # calculating sin and cosin of the aspect layer
    exp = "$out = cos($a)"
    grass.mapcalc(exp, out="coseno", a="aspectslutmp", quiet=True)
    exp = "$out = sin($a)"
    grass.mapcalc(exp, out="seno", a="aspectslutmp", quiet=True)

    grass.run_command(
        "g.remove", type="raster", name="aspectslutmp", flags="f", quiet=True
    )
    rm_rasters.append("coseno")
    rm_rasters.append("seno")

    # setting counters and base layers for the next "while loop"
    counter = 0
    last_counter = counter - 1
    control = 1
    i = grass.raster_info(dem)
    # control_lastrun = int(i["cells"])
    grass.mapcalc("$out = null()", out="slu_r_0", quiet=True)
    grass.mapcalc("$out = null()", out="cvar_0", quiet=True)
    grass.mapcalc("$out = null()", out="count_0", quiet=True)
    grass.mapcalc("$out = 1", out="slu_r_todo", quiet=True)
    grass.mapcalc("$out = 1", out="slu_r_todo_0", quiet=True)
    rm_rasters.append("slu_r_0")
    rm_rasters.append("cvar_0")
    rm_rasters.append("count_0")
    rm_rasters.append("slu_r_todo")
    rm_rasters.append("slu_r_todo_0")

    # starting the loop. The loop stops when: there are no halfbasins extracted (control=0),
    # OR the number of allowed iteration is exceeded (counter>maxiter) OR the threshold
    # (cells) is greather/equal than the reduction factor (thc >= red) otherwise int(thc-thc/red)
    # remains equal to thc
    while control > 0 and counter < maxiter and thc >= red:

        # generating half-basins
        # grass.run_command('r.watershed', elevation=dem, hbasin='slu_r_tmp', thresh=thc, flags='abs', overwrite=True, quiet=True)
        grass.run_command(
            "g.remove", type="raster", name="slu_r_tmp", flags="f", quiet=True
        )
        grass.run_command(
            "r.watershed",
            elevation=dem,
            hbasin="slu_r_tmp",
            thresh=thc,
            flags="ab",
            quiet=True,
        )
        rm_rasters.append("slu_r_tmp")

        # remove output map from previous run
        grass.run_command(
            "g.remove", type="raster", name="slu_r", flags="f", quiet=True
        )
        if options["plainsmap"]:
            exp = "$out = if(isnull($a), $b, null())"
            grass.mapcalc(exp, out="slu_r", a=plains, b="slu_r_tmp", quiet=True)
        else:
            grass.run_command("g.copy", rast="slu_r_tmp,slu_r", quiet=True)

        rm_rasters.append("slu_r")

        grass.run_command("r.mask", flags="r", quiet=True)
        exp = "$out = if(isnull($mask), null(), 1)"
        grass.mapcalc(exp, out="MASK", mask="slu_r_todo", quiet=True)
        grass.run_command(
            "r.stats.zonal",
            base="slu_r",
            cover="coseno",
            method="count",
            output="count",
            quiet=True,
        )
        rm_rasters.append("count")
        grass.run_command(
            "r.stats.zonal",
            base="slu_r",
            cover="coseno",
            method="sum",
            output="sumcoseno",
            quiet=True,
        )
        rm_rasters.append("sumcoseno")
        grass.run_command(
            "r.stats.zonal",
            base="slu_r",
            cover="seno",
            method="sum",
            output="sumseno",
            quiet=True,
        )
        rm_rasters.append("sumseno")

        # creating, for each half-basin, the layer where the circular variance is stored (cvar).
        # Circular variance is 1-R/n, R is the magnitude of the vectorial sum of all the unit
        # vectors of the aspect layer in each polygon and n is the number of unit vectors (and
        # cells) involved in the sum
        exp = "$out = 1 - ((sqrt(($a)^2 + ($b)^2)) / $c)"
        grass.mapcalc(
            exp,
            out="cvar",
            a="sumseno",
            b="sumcoseno",
            c="count",
            quiet=True,
        )
        rm_rasters.append("cvar")

        grass.run_command("r.mask", flags="r", quiet=True)

        # selecting half-basins where area is larger than the minimum area and the average
        # unit vector is smaller than the unit vector threshold
        grass.run_command(
            "g.remove", type="raster", name="slu_r_todo", flags="f", quiet=True
        )
        if options["areamax"]:
            exp = "$out = if($a > $f || ($a > $b && $c > $d), $g, null())"
            grass.mapcalc(
                exp,
                out="slu_r_todo",
                a="count",
                b=aminc,
                c="cvar",
                d=cvarmin,
                g="slu_r",
                f=amaxc,
                quiet=True,
            )
            grass.run_command(
                "g.copy",
                rast=("count,count_prova_%d") % counter,
                quiet=True,
            )
            rm_rasters.append("count_prova_%d" % counter)
            # exp = "$out = if($a>$b,1,null())"
            # grass.mapcalc(exp, out = "slu_r_large"+str(counter), a = "count", b = amaxc , overwrite=True, quiet=True)
        else:
            exp = "$out = if($a > $b && $c > $d, $g, null())"
            grass.mapcalc(
                exp,
                out="slu_r_todo",
                a="count",
                b=aminc,
                c="cvar",
                d=cvarmin,
                g="slu_r",
                quiet=True,
            )

        # checking that there actually are half-basins with area greater than areamin
        # and circular variance greater than cvarmin. otherwise the loop exits
        s = grass.read_command("r.univar", flags="g", map="slu_r_todo", quiet=True)
        kv = grass.parse_key_val(s)
        # ivan
        # if there are any non-NULL cells
        if int(kv["n"]) > 0:
            # if kv.has_key("n"):
            # increasing counter
            last_counter = counter
            counter = counter + 1
            # patching the current half-basins, cvar and counter that were not selected
            # in the previous steps with those that come from the previous step of the loop
            grass.run_command(
                "g.copy", rast=("slu_r_todo,slu_r_todo_%d") % counter, quiet=True
            )
            rm_rasters.append("slu_r_todo_%s" % counter)
            grass.run_command("r.mask", raster="slu_r_todo", flags="i", quiet=True)
            grass.run_command(
                "r.patch",
                input=("slu_r_" + str(last_counter), "slu_r"),
                output="slu_r_" + str(counter),
                quiet=True,
            )
            rm_rasters.append("slu_r_%d" % counter)
            grass.run_command(
                "g.copy",
                rast=("slu_r_" + str(counter), "slu_r_prova_" + str(counter)),
                quiet=True,
            )
            rm_rasters.append("slu_r_prova_%d" % counter)
            grass.run_command(
                "r.patch",
                input=("cvar_" + str(last_counter), "cvar"),
                output="cvar_" + str(counter),
                quiet=True,
            )
            rm_rasters.append("cvar_%d" % counter)
            grass.run_command(
                "r.patch",
                input=("count_" + str(last_counter), "count"),
                output="count_" + str(counter),
                quiet=True,
            )
            rm_rasters.append("count_%d" % counter)

            grass.run_command("r.mask", flags="r", quiet=True)

            # rejecting partition if average area of new half-basins is less than amin;
            # not effective on large areas, if areamax is present
            if counter > 0:  # counter is always > 0 ?
                if options["areamax"]:
                    # this block does not make sense,
                    # count_prova_<last_counter> has been calculated above
                    if counter == 1:
                        grass.mapcalc("$out = 1", out="count_prova_0", quiet=True)
                    exp = "$out = if($a > $b, 1, null())"
                    # this
                    # a="count_prova_" + str(last_counter - 1)
                    # can not work if counter = 1 and last_counter = 0
                    # because the result would be "count_prova_-1"
                    # changed to
                    # a="count_prova_" + str(last_counter)
                    grass.mapcalc(
                        exp,
                        out="slu_r_large" + str(counter),
                        a="count_prova_" + str(last_counter),
                        b=amaxc,
                        quiet=True,
                    )
                    rm_rasters.append("slu_r_large%d" % counter)
                    exp = "$out = if(isnull($a), $b, null())"
                    grass.mapcalc(
                        exp,
                        out="MASK",
                        a="slu_r_large" + str(counter),
                        b="slu_r_" + str(counter),
                        quiet=True,
                    )
                    grass.run_command(
                        "g.copy",
                        rast="MASK,mask" + str(counter),
                        quiet=True,
                    )
                else:
                    grass.run_command(
                        "r.mask", raster="slu_r_" + str(counter), quiet=True
                    )

                z = grass.read_command(
                    "r.univar",
                    flags="g",
                    map="slu_r_todo_" + str(last_counter),
                    quiet=True,
                )
                kvz = grass.parse_key_val(z)

                #                grass.message(("Univar: %s") % kvz )
                # ivan
                # if kvz.has_key("n"):
                # if there are any non-NULL cells
                if int(kvz["n"]) > 0:
                    en = int(kvz["n"])
                    grass.message(("Univar: %s") % en)
                    if en > 0:
                        grass.run_command(
                            "r.statistics",
                            base="slu_r_todo_" + str(last_counter),
                            cover="slu_r",
                            method="diversity",
                            output="slu_diversity_" + str(counter),
                            quiet=True,
                        )
                        rm_rasters.append("slu_diversity_%d" % counter)
                        #                        grass.run_command('r.univar', map='slu_r')
                        #                        grass.run_command('r.univar', map='slu_r_todo_'+str(last_counter))
                        #                        grass.run_command('r.univar', map='slu_diversity_'+str(counter))
                        grass.run_command(
                            "r.stats.zonal",
                            base="slu_r_todo_" + str(last_counter),
                            cover="coseno",
                            method="count",
                            output="todocount_" + str(counter),
                            quiet=True,
                        )
                        rm_rasters.append("todocount_%d" % counter)
                        exp = "$out = $d / $a"
                        grass.mapcalc(
                            exp,
                            out="slu_r_test_" + str(counter),
                            a="slu_diversity_" + str(counter),
                            d="todocount_" + str(counter),
                            quiet=True,
                        )
                        rm_rasters.append("slu_r_test_%d" % counter)
                        exp = "$out = if($d < $e, $c, null())"
                        grass.mapcalc(
                            exp,
                            out="slu_r_corr_" + str(counter),
                            b="slu_r_" + str(counter),
                            c="slu_r_todo_" + str(last_counter),
                            d="slu_r_test_" + str(counter),
                            e=aminc,
                            quiet=True,
                        )
                        rm_rasters.append("slu_r_corr_%d" % counter)
                        grass.run_command("r.mask", flags="r", quiet=True)

                        grass.run_command(
                            "g.rename",
                            rast=f"slu_r_{counter},slu_r_tmp_{counter}",
                            quiet=True,
                        )
                        rm_rasters.append("slu_r_tmp_%d" % counter)
                        grass.run_command(
                            "r.patch",
                            input=(
                                "slu_r_corr_" + str(counter),
                                "slu_r_tmp_" + str(counter),
                            ),
                            output="slu_r_" + str(counter),
                            quiet=True,
                        )
                        grass.run_command(
                            "g.remove",
                            type="raster",
                            name=f"slu_r_tmp_{counter}",
                            flags="f",
                            quiet=True,
                        )
                        rm_rasters.append("slu_r_%d" % counter)
                    else:
                        grass.run_command("r.mask", flags="r", quiet=True)
                else:
                    grass.run_command("r.mask", flags="r", quiet=True)

            control = int(kv["n"])
            thc = int(thc - thc / red)
            thhect = thc * nsres * ewres / 10000
            grass.message(("Threshold (hectars) is: %s") % thhect)
            grass.message(
                ("No. of cells to be still classified as SLU is: %s. Loop done: %s")
                % (control, counter)
            )
        else:
            # exit the loop
            grass.message(("Nothing to do, ready to write the outputs"))
            control = 0

    # depending on how the while loop is exited the slu_r_$counter may have some small holes. Here we fill them.
    exp = "$out = if(isnull($a), $b, $a)"
    grass.mapcalc(
        exp,
        out="slu_r_final",
        a="slu_r_" + str(counter),
        b="slu_r",
        quiet=True,
    )
    exp = "$out = $a"
    grass.mapcalc(
        exp,
        out="cvar_final",
        a="cvar_" + str(last_counter),
        quiet=True,
    )
    exp = "$out = $a"
    grass.mapcalc(
        exp,
        out="count_final",
        a="count_" + str(last_counter),
        quiet=True,
    )

    # preparing the outputs
    exp = "$out = $a"
    grass.mapcalc(exp, out="slumap_1", a="slu_r_final", quiet=True)
    # add areas where DEM exists, and SUs do not exist
    if options["plainsmap"]:
        exp = "$out = if(isnull($a), if(isnull($b), if(isnull($c), null(), 1), null()), $a)"
        grass.mapcalc(
            exp,
            a="slumap_1",
            b=plains,
            c=dem,
            out="slumap_2",
            quiet=True,
        )
    else:
        exp = "$out = if(isnull($a), if(isnull($c), null(), 1), $a)"
        grass.mapcalc(exp, a="slumap_1", c=dem, out="slumap_2", quiet=True)
    grass.run_command("r.clump", input="slumap_2", output=slumap, quiet=True)
    grass.run_command("g.remove", type="raster", name="slumap_1", flags="f", quiet=True)
    grass.run_command("g.remove", type="raster", name="slumap_2", flags="f", quiet=True)

    grass.run_command("r.colors", map="slu_r_final", color="random", quiet=True)

    if circvarmap:
        grass.message("Writing out %s" % circvarmap)
        exp = "$out = $a"
        grass.mapcalc(exp, out=circvarmap, a="cvar_final", quiet=True)
    if areamap:
        grass.message("Writing out %s" % areamap)
        exp = "$out = $a * $b * $c"
        grass.mapcalc(
            exp,
            out=areamap,
            a="count_final",
            b=nsres,
            c=ewres,
            quiet=True,
        )

    # clean up
    grass.run_command("r.mask", flags="r", quiet=True)
    grass.run_command(
        "g.remove", type="raster", name="seno,coseno", flags="f", quiet=True
    )

    grass.run_command(
        "g.remove", type="raster", name="slu_r_todo", flags="f", quiet=True
    )
    grass.run_command(
        "g.remove",
        type="raster",
        name="slu_r_final,cvar_final,count_final,slu_r_todo_final",
        flags="f",
        quiet=True,
    )

    for i in range(counter + 1):
        mapname = "slu_r_" + str(i)
        grass.run_command(
            "g.remove", type="raster", name=mapname, flags="f", quiet=True
        )
        mapname = "cvar_" + str(i)
        grass.run_command(
            "g.remove", type="raster", name=mapname, flags="f", quiet=True
        )
        mapname = "count_" + str(i)
        grass.run_command(
            "g.remove", type="raster", name=mapname, flags="f", quiet=True
        )
        mapname = "slu_r_todo_" + str(i)


def clean_method_3(input_vect, output_vect, minarea):
    region = grass.region()
    nsres = region["nsres"]
    ewres = region["ewres"]
    smarea = 10 * nsres * ewres

    grass.run_command(
        "v.clean",
        input=input_vect,
        output="slu_clean",
        tool="rmarea",
        threshold=smarea,
        quiet=True,
    )
    grass.run_command(
        "v.db.addcolumn",
        map="slu_clean",
        columns="area integer, perimetro integer",
        quiet=True,
    )
    grass.run_command(
        "v.to.db", map="slu_clean", option="area", columns="area", quiet=True
    )
    grass.run_command(
        "v.to.db", map="slu_clean", option="perimeter", columns="perimetro", quiet=True
    )
    grass.run_command(
        "v.db.droprow",
        input="slu_clean",
        where="area is null",
        output="slu_area",
        quiet=True,
    )
    rm_vectors.append("slu_clean")

    key = grass.vector_db(map="slu_area")[1]["key"]
    lista = grass.read_command(
        "v.db.select", map="slu_area", columns=key, where=f"area <= {minarea}", flags="c"
    )
    lista = lista.splitlines()
    buchi = ",".join(lista)
    totalebuchi = len(lista)
    grass.run_command(
        "v.extract",
        input="slu_area",
        cats=buchi,
        output="slu_buchi",
        type="area",
        quiet=True,
    )
    grass.run_command(
        "v.to.rast", input="slu_buchi", output="slu_buchi", use="cat" - -o - -q
    )

    key = grass.vector_db(map="slu_area")[1]["key"]
    lista = grass.read_command(
        "v.db.select", map="slu_area", columns=key, where=f"area > {minarea}", flags="c"
    )
    lista = lista.splitlines()
    nobuchi = ",".join(lista)
    grass.run_command(
        "v.extract",
        input="slu_area",
        cats=nobuchi,
        output="slu_nobuchi",
        type="area",
        quiet=True,
    )
    grass.run_command(
        "v.to.rast", input="slu_nobuchi", output="slu_nobuchi", use="cat" - -o - -q
    )
    grass.run_command(
        "v.category",
        input="slu_area",
        output="slu_bordi",
        layer=2,
        type="boundary",
        option="add",
        quiet=True,
    )
    grass.run_command(
        "v.db.addtable",
        map="slu_bordi",
        layer=2,
        columns="left integer,right integer,lunghezza integer",
        quiet=True,
    )
    grass.run_command(
        "v.to.db",
        map="slu_bordi",
        option="sides",
        columns="left,right",
        layer=2,
        type="boundary",
        quiet=True,
    )
    grass.run_command(
        "v.to.db",
        map="slu_bordi",
        option="length",
        columns="lunghezza",
        layer=2,
        type="boundary",
        quiet=True,
    )
    grass.run_command(
        "v.to.rast", input="slu_area", output="slu_area", use="cat", quiet=True
    )
    # TODO: different names than coseno and seno, these rasters are already created
    grass.mapcalc("coseno = cos(aspect_slu)", quiet=True)
    grass.mapcalc("seno = sin(aspect_slu)", quiet=True)

    grass.run_command(
        "r.stats.zonal",
        base="slu_area",
        cover="coseno",
        method="count",
        output="count",
        quiet=True,
    )
    grass.run_command(
        "r.stats.zonal",
        base="slu_area",
        cover="coseno",
        method="sum",
        output="cumcos",
        quiet=True,
    )
    grass.run_command(
        "r.stats.zonal",
        base="slu_area",
        cover="seno",
        method="sum",
        output="sumsin",
        quiet=True,
    )

    grass.mapcalc("cos_medio = sumcos / count", quiet=True)
    grass.mapcalc("sin_medio = sumsin / count", quiet=True)
    grass.run_command(
        "r.to.vect", input="cos_medio", output="cos_medio", type="area", quiet=True
    )
    grass.run_command(
        "r.to.vect", input="sin_medio", output="sin_medio", type="area", quiet=True
    )
    grass.run_command(
        "v.overlay",
        ainput="slu_area",
        binput="cos_medio",
        operator="and",
        atype="area",
        btype="area",
        output="cos_medio_over",
        quiet=True,
    )
    grass.run_command(
        "v.overlay",
        ainput="slu_area",
        binput="sin_medio",
        operator="and",
        atype="area",
        btype="area",
        output="sin_medio_over",
        quiet=True,
    )

    pulire = grass.read_command(
        "v.category", input="slu_buchi", option="print", quiet=True
    )
    pulire = pulire.splitlines()

    grass.run_command("g.copy", vector=f"slu_area,{output_vect}", quiet=True)

    ico = 1
    for i in pulire:
        inti = int(i)
        lista1 = grass.read_command(
            "db.select",
            sql=f"select b2.right from slu_bordi_2 b2 where b2.left = {i} and b2.right <> -1",
            flags="c",
        )
        lista2 = grass.read_command(
            "db.select",
            sql=f"select b2.left from slu_bordi_2 b2 where b2.left <> -1 and b2.right = {inti}",
            flags="c",
        )
        vicini = lista1.splitlines()
        vicini.extend(lista2.splitlines)
        vicini = sorted(set(vicini))
        if len(vicini) > 0:
            grass.message(
                f" --- --- -- buco numero {ico} di {totalebuchi}, cat: {i}, vicini: {vicini}"
            )
            ico = ico + 1
            grass.message(f"vicini: {vicini}")
            grass.run_command(
                "v.extract",
                input=output_vect,
                cats=",".join(vicini),
                output="intorno",
                type="area",
                quiet=True,
            )
            chk_intorno = grass.read_command(
                "v.category",
                input="intorno",
                type="centroid",
                option="print",
                quiet=True,
            )
            chk_intorno = chk_intorno.splitlines()
            if len(chk_intorno) > 0:
                # potrei voler cambiare questo perche' quando ci sono buchi contigui fa un po' di casino
                grass.run_command(
                    "v.overlay",
                    ainput="intorno",
                    binput="slu_nobuchi",
                    output="intorno_OK",
                    atype="area",
                    btype="area",
                    olayer="0,1,0",
                    operator="and",
                    quiet=True,
                )

                cos_buco = grass.read_command(
                    "v.db.select",
                    map="cos_medio_over",
                    where=f"a_cat={i}",
                    columns="b_value",
                    flags="c",
                    quiet=True,
                )
                sin_buco = grass.read_command(
                    "v.db.select",
                    map="sin_medio_over",
                    where=f"a_cat={i}",
                    columns="b_value",
                    flags="c",
                    quiet=True,
                )
                grass.message(f"buco cat {i}: cos={cos_buco} sin={sin_buco}")
                massimo = -10000
                jmax = 0
                loop = grass.read_command(
                    "v.category", input="intorno_OK", option="print", quiet=True
                )
                loop = loop.splitlines()
                for j in loop:
                    j = int(j)
                    cos_j = grass.read_command(
                        "v.db.select",
                        map="cos_medio_over",
                        where=f"a_cat={j}",
                        columns="b_value",
                        flags="c",
                        quiet=True,
                    )
                    sin_j = grass.read_command(
                        "v.db.select",
                        map="sin_medio_over",
                        where=f"a_cat={j}",
                        columns="b_value",
                        flags="c",
                        quiet=True,
                    )
                    dotpr = (
                        float(cos_buco) * float(cos_j) + float(sin_buco) * float(sin_j)
                    ) * 10000
                    if dotpr >= massimo and dotpr > 0:
                        massimo = dotpr
                        jmax = j

                    grass.message(
                        f"i: {i} j: {j} cos_j: {cos_j} sin_j: {sin_j} dotpr: {dotpr} jmax: {jmax}"
                    )

                grass.message(f"massimo: {massimo} per j={jmax}")
                if jmax > 0:
                    lunghezza = grass.read_command(
                        "db.select",
                        sql=f"select b2.lunghezza from slu_bordi_2 b2 where (b2.left={i} and b2.right={jmax}) or (b2.left={jmax} and b2.right={i})",
                        flags="c",
                        quiet=True,
                    )
                    perimetro = grass.read_command(
                        "v.db.select",
                        map="slu_clean",
                        columns="perimetro",
                        where=f"cat={i}",
                        flags="c",
                        quiet=True,
                    )
                    lunghezza = float(lunghezza)
                    perimetro = float(perimetro)
                    if lunghezza > 0 and perimetro > 0:
                        frazione = lunghezza / perimetro * 10000
                        if frazione > 500:
                            grass.message(
                                f"lungh: {lunghezza}; perim: {perimetro}; fract: {frazione}"
                            )
                            grass.run_command(
                                "v.extract",
                                input=output_vect,
                                output="slu_i",
                                cats=f"{i},{jmax}",
                                new={jmax},
                                flags="d",
                                quiet=True,
                            )
                            grass.run_command(
                                "v.overlay",
                                ainput=output_vect,
                                binput="slu_i",
                                atype="area",
                                btype="area",
                                operator="not",
                                olayer="0,1,0",
                                output="slu_j",
                                quiet=True,
                            )
                            grass.run_command(
                                "v.overlay",
                                ainput="slu_i",
                                binput="slu_j",
                                atype="area",
                                btype="area",
                                operator="or",
                                output="slu_k",
                                olayer="1,0,0",
                                quiet=True,
                            )
                            grass.run_command(
                                "v.db.addcolumn",
                                map="slu_k",
                                column="newcat integer",
                                quiet=True,
                            )
                            grass.run_command(
                                "v.db.update",
                                map="slu_k",
                                layer=1,
                                column="newcat",
                                qcolumn="a_cat",
                                where="a_cat is not null",
                                quiet=True,
                            )
                            grass.run_command(
                                "v.db.update",
                                map="slu_k",
                                layer=1,
                                column="newcat",
                                qcolumn="b_cat",
                                where="b_cat is not null",
                                quiet=True,
                            )
                            grass.run_command(
                                "v.reclass",
                                input="slu_k",
                                output=output_vect,
                                column="newcat",
                                quiet=True,
                            )
                            grass.run_command(
                                "v.db.addtable", map=output_vect, quiet=True
                            )
                            grass.run_command(
                                "v.db.addcolumn",
                                map=output_vect,
                                columns="area integer",
                                quiet=True,
                            )
                            grass.run_command(
                                "v.to.db",
                                map=output_vect,
                                option="area",
                                columns="area",
                                quiet=True,
                            )
                            grass.run_command(
                                "g.remove",
                                type="vector",
                                name="slu_i,slu_j,slu_k",
                                flags="f",
                                quiet=True,
                            )
                    # lunghezza and perimetro
                # jmax
                grass.run_command(
                    "g.remove", type="vector", name="intorno_OK", flags="f", quiet=True
                )
            # chk_category
            grass.run_command(
                "g.remove", type="vector", name="intorno", flags="f", quiet=True
            )
        # vicini
    # pulire


def clean_small_areas(dem, slumap, plains, cleansize=-1, slumapclean=None):
    if not flags["n"]:
        if not flags["m"]:
            grass.message(
                " -- we want QUICK cleaning of small-sized areas: METHOD 1 --"
            )

        exp = "$out = if(isnull($mask), null(), 1)"
        grass.mapcalc(exp, out="MASK", mask=dem, quiet=True)
        areamap = "areamap"

        grass.run_command("r.clump", input=slumap, output="slu_clump", quiet=True)
        rm_rasters.append("slu_clump")

        grass.run_command(
            "r.stats.zonal",
            base="slu_clump",
            cover="slu_clump",
            method="count",
            output="slu_count",
            quiet=True,
        )
        rm_rasters.append("slu_count")

        exp = "$out = $a * $b * $c"
        grass.mapcalc(
            exp,
            out=areamap,
            a="slu_count",
            b=nsres,
            c=ewres,
            quiet=True,
        )
        rm_rasters.append(areamap)

        exp = "$out = if($a > $b, $c, null())"
        grass.mapcalc(
            exp,
            out="slu_r_clean",
            a=areamap,
            b=cleansize,
            c=slumap,
            quiet=True,
        )
        rm_rasters.append("slu_r_clean")

        cleansize = cleansize / (nsres * ewres)
        growdist = int((10 * cleansize / 3.14) ** 0.5)
        grass.run_command(
            "r.grow",
            input="slu_r_clean",
            output="slu_r_grow",
            radius=growdist,
            quiet=True,
        )
        rm_rasters.append("slu_r_grow")

    if flags["m"]:
        grass.message(" -- we want QUICK cleaning of small-sized areas: METHOD 2 --")
        input = "slu_r_grow"
        output = "slu_no_stripes"
        grass.run_command(
            "r.neighbors",
            input=input,
            output="slu_diversity",
            method="diversity",
            size=5,
            quiet=True,
        )
        rm_rasters.append("slu_diversity")

        exp = "$out = if($a == 1, 1, null())"
        grass.mapcalc(
            exp,
            out="slu_diversity_nobordi",
            a="slu_diversity",
            quiet=True,
        )
        rm_rasters.append("slu_diversity_nobordi")

        grass.run_command(
            "r.grow",
            input="slu_diversity_nobordi",
            output="slu_diversity_nobordi_grow",
            radius=1.01,
            quiet=True,
        )
        rm_rasters.append("slu_diversity_nobordi_grow")

        exp = "$out = if(isnull($a), null(), $b)"
        grass.mapcalc(
            exp,
            out="slu_finale_nobordi",
            a="slu_diversity_nobordi_grow",
            b=input,
            quiet=True,
        )
        rm_rasters.append("slu_finale_nobordi")

        # TODO: fixed radius ?
        grass.run_command(
            "r.grow",
            input="slu_finale_nobordi",
            output=output,
            radius=1000,
            quiet=True,
        )
        rm_rasters.append(output)

        exp = "$out = int($a)"
        # TODO: input="slu_r_grow" exists already
        grass.mapcalc(exp, out=input, a=output, quiet=True)
        grass.run_command(
            "g.remove",
            type="raster",
            name="slu_diversity,slu_diversity_nobordi,slu_diversity_nobordi_grow,slu_finale_nobordi,slu_no_stripes",
            flags="f",
            quiet=True,
        )

    if flags["n"]:
        grass.message(" -- we want DETAILED cleaning of small-sized areas: METHOD 3 --")
        grass.run_command(
            "r.to.vect",
            input=slumap,
            output="slu_v_grow",
            type="area",
            quiet=True,
        )
        rm_vectors.append("slu_v_grow")
        # TODO: include in new fn clean_method_3
        grass.run_command(
            "r.slope.aspect",
            elevation=dem,
            aspect="aspect_slu",
            quiet=True,
        )
        rm_rasters.append("aspect_slu")
        # os.system(
        #    "/home/alvioli/dem_globo_srtm/range/script04/slu_code/clean_method_3.sh slu_v_grow vect2 %s"
        #    % cleansize
        # )
        clean_method_3("slu_v_grow", "vect2", cleansize)

        # applying method 2 at the end
        grass.run_command(
            "v.to.rast",
            input="vect2",
            output="rast2",
            use="cat",
            quiet=True,
        )
        rm_rasters.append("rast2")

        input = "rast2"
        output = "slu_r_grow"
        grass.run_command(
            "r.neighbors",
            input=input,
            output="slu_diversity",
            method="diversity",
            size=5,
            quiet=True,
        )
        rm_rasters.append("slu_diversity")

        exp = "$out = if($a == 1, 1, null())"
        grass.mapcalc(
            exp,
            out="slu_diversity_nobordi",
            a="slu_diversity",
            quiet=True,
        )
        rm_rasters.append("slu_diversity_nobordi")

        grass.run_command(
            "r.grow",
            input="slu_diversity_nobordi",
            output="slu_diversity_nobordi_grow",
            radius=1.01,
            quiet=True,
        )
        rm_rasters.append("slu_diversity_nobordi_grow")

        exp = "$out = if(isnull($a), null(), $b)"
        grass.mapcalc(
            exp,
            out="slu_finale_nobordi",
            a="slu_diversity_nobordi_grow",
            b=input,
            quiet=True,
        )
        rm_rasters.append("slu_finale_nobordi")

        grass.run_command(
            "r.grow",
            input="slu_finale_nobordi",
            output=output,
            radius=1000,
            quiet=True,
        )
        rm_rasters.append(output)

        exp = "$out = int($a)"
        grass.mapcalc(exp, out=input, a=output, quiet=True)
        rm_rasters.append(input)

        grass.run_command(
            "g.remove",
            type="raster",
            name="slu_diversity,slu_diversity_nobordi,slu_diversity_nobordi_grow,slu_finale_nobordi,rast2",
            flags="f",
            quiet=True,
        )

        grass.run_command(
            "g.remove",
            type="vector",
            name="slu_v_grow,vect2",
            flags="f",
            quiet=True,
        )

    if options["plainsmap"]:
        exp = "$out = if(isnull($b), if(isnull($c), null(), int($a)), null())"
        grass.mapcalc(
            exp,
            out=slumapclean,
            a="slu_r_grow",
            b=plains,
            c=dem,
            quiet=True,
        )
    else:
        exp = "$out = if(isnull($c), null(), int($a))"
        grass.mapcalc(exp, out=slumapclean, a="slu_r_grow", c=dem, quiet=True)
    grass.run_command("r.colors", map=slumapclean, color="random", quiet=True)


def main():
    global rm_rasters, rm_vectors
    global thc

    dem = options["demmap"]
    plains = None
    if options["plainsmap"]:
        plains = options["plainsmap"]
    slumap = options["slumap"]
    circvarmap = options["circvarmap"]
    areamap = options["areamap"]
    th = float(options["thresh"])
    amin = float(options["areamin"])
    cvarmin = float(options["cvmin"])
    red = int(options["rf"])
    maxiter = int(options["maxiteration"])

    cleansize = -1
    slumapclean = None
    if options["cleansize"]:
        cleansize = int(options["cleansize"])
        if options["slumapclean"]:
            slumapclean = options["slumapclean"]
        else:
            grass.fatal("When cleansize is provided, slumapclean is mandatory.")
        if flags["m"] and flags["n"]:
            grass.fatal(
                "When cleansize is provided, only one between m and n can be specified."
            )

    slope_units(
        dem, plains, slumap, circvarmap, areamap, th, amin, cvarmin, red, maxiter
    )

    if cleansize > 0:
        clean_small_areas(dem, slumap, plains, cleansize, slumapclean)


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
