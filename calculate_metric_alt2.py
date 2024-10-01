#!/usr/bin/env python3


import grass.script as grass


# HINT: "echo" of exact GRASS GIS command was not translated
# HINT: quiet / superquiet flags were not translated


def calculate_metric_alt2(
    basin, dem, sucl, thr_clean, uid, areamin, cvmin, res, outfile
):
    """Calculate metric"""

    grass.run_command("g.region", vect=basin, align=dem)
    grass.run_command("r.mask", flags="r")
    grass.run_command("r.mask", vect=basin, overwrite=True)
    grass.run_command("g.remove", type="vector", name="su_segm", flags="f")
    grass.run_command("g.copy", vect=f"{sucl},su_segm", overwrite=True)
    grass.run_command(
        "v.db.dropcolumn", map="su_segm", columns="value,label,area"
    )
    grass.run_command("v.db.addcolumn", map="su_segm", columns="area real")
    grass.run_command(
        "v.to.db", map="su_segm", option="area", columns="area", overwrite=True
    )
    grass.run_command(
        "db.execute", sql=f"delete from su_segm where area<{thr_clean}"
    )
    grass.run_command(
        "v.clean",
        input="su_segm",
        type="area",
        output="su_ok",
        tool="rmarea",
        threshold=thr_clean,
        overwrite=True,
    )
    grass.run_command("g.rename", vect="su_ok,su_segm", overwrite=True)

    ##############################################################
    #
    # Calculation of the procedure proposed in
    # [1] Espindola et al., Int. J. Remote Sens. 27, 3035-3040 (2006)
    # and adapted for our case in
    # [2] Alvioli, Marchesini et al., Geosci. Mod. Dev. 9, 3975-3991 (2016)
    # and
    # [3] Alvioli, Guzzetti, Marchesini, Geomorphology 358, 107124 (2020)
    # [4] Alvioli, Marchesini et al., Journal of Maps, 18, 300-313 (2021)
    #

    grass.run_command("r.mask", flags="r")
    grass.run_command("r.mask", vect="su_segm", overwrite=True)
    grass.run_command(
        "v.to.rast",
        input="su_segm",
        output="su_segm",
        use="cat",
        overwrite=True,
    )

    onecell = res * res

    grass.message(_("r.slope.aspect ..."))
    grass.run_command(
        "r.slope.aspect", elevation=dem, aspect="aspect", overwrite=True
    )
    grass.run_command("r.colors", map="aspect", color="aspectcolr")
    grass.message(_("circular variance ..."))
    # aspect in degrees - mapcalc's cos(x) wants x in degrees -- OK
    grass.run_command(
        "r.mapcalc", expression="cos=cos(aspect)", overwrite=True
    )
    # aspect in degrees - mapcalc's sin(x) wants x in degrees -- OK
    grass.run_command(
        "r.mapcalc", expression="sin=sin(aspect)", overwrite=True
    )

    grass.message(_("r.statistics 1 ..."))
    grass.run_command(
        "r.stats.zonal",
        base="su_segm",
        cover="su_segm",
        method="count",
        output="numero",
        overwrite=True,
    )
    grass.message(_("r.statistics 2 ..."))
    grass.run_command(
        "r.stats.zonal",
        base="su_segm",
        cover="cos",
        method="sum",
        output="sumcos",
        overwrite=True,
    )
    grass.message(_("r.statistics 3 ..."))
    grass.run_command(
        "r.stats.zonal",
        base="su_segm",
        cover="sin",
        method="sum",
        output="sumsin",
        overwrite=True,
    )
    grass.run_command(
        "r.mapcalc",
        expression="v_i = 1-((sqrt((sumsin)^2 + (sumcos)^2))/numero)",
        overwrite=True,
    )

    grass.message(_("r.statistics 4 ..."))
    grass.run_command(
        "r.mapcalc", expression="base = int($dem/$dem)", overwrite=True
    )
    grass.run_command(
        "r.stats.zonal",
        base="base",
        cover="cos",
        method="sum",
        output="sumcos_all",
        overwrite=True,
    )
    grass.run_command(
        "r.stats.zonal",
        base="base",
        cover="sin",
        method="sum",
        output="sumsin_all",
        overwrite=True,
    )

    grass.message(_("v.category edges ..."))
    grass.run_command(
        "g.remove", type="vector", name="su_segm_edges", flags="f"
    )
    grass.run_command(
        "v.category",
        input="su_segm",
        output="su_segm_edges",
        layer=2,
        type="boundary",
        option="add",
        overwrite=True,
    )
    grass.run_command(
        "v.db.addtable",
        map="su_segm_edges",
        layer=2,
        col="left integer,right integer,lunghezza real",
    )
    grass.run_command(
        "v.to.db",
        map="su_segm_edges",
        option="sides",
        columns="left,right",
        layer=2,
        type="boundary",
        overwrite=True,
    )
    grass.run_command(
        "v.to.db",
        map="su_segm_edges",
        option="length",
        columns="lunghezza",
        layer=2,
        type="boundary",
        overwrite=True,
    )

    # mapcalc's atan(x,y) returns result in degrees ..
    grass.run_command(
        "r.mapcalc", expression="a_i = atan(sumsin,sumcos)", overwrite=True
    )
    grass.run_command(
        "r.mapcalc",
        expression="a_all = atan(sumsin_all,sumcos_all)",
        overwrite=True,
    )
    # .. but we need maps in radians, as needed by bc in bash
    pii = 3.14159265359
    grass.run_command(
        "r.mapcalc", expression=f"a_irad = a_i*{pii}/180", overwrite=True
    )
    grass.run_command(
        "r.mapcalc", expression=f"a_allrad = a_all*{pii}/180", overwrite=True
    )

    grass.message(_("v.rast.stats 1 ..."))
    grass.run_command(
        "v.rast.stats",
        map="su_segm_edges",
        raster="a_irad",
        column_prefix="a_i",
        method="average",
    )
    grass.run_command(
        "v.db.renamecolumn", map="su_segm_edges", column="a_i_average,a_i"
    )

    grass.message(_("v.rast.stats 2 ..."))
    grass.run_command(
        "v.rast.stats",
        map="su_segm_edges",
        raster="a_allrad",
        column_prefix="a_all",
        method="average",
    )
    grass.run_command(
        "v.db.renamecolumn", map="su_segm_edges", column="a_all_average,a_all"
    )

    grass.message(_("v.rast.stats 3 ..."))
    grass.run_command(
        "v.rast.stats",
        map="su_segm_edges",
        raster="v_i",
        column_prefix="v_i",
        method="average",
    )
    grass.run_command(
        "v.db.renamecolumn", map="su_segm_edges", column="v_i_average,v_i"
    )

    # fix problems
    grass.run_command(
        "db.execute",
        sql=(
            "update su_segm_edges set v_i = (select coalesce(0,v_i) "
            "as v_i from su_segm_edges) where v_i is null"
        ),
    )
    grass.run_command(
        "v.extract",
        input="su_segm_edges",
        where="v_i>=0",
        output="su_segm",
        overwrite=True,
    )
    grass.run_command(
        "v.db.addtable", map="su_segm", table="su_segm_edges_2", layer=2
    )

    # calcolo V - metodo RASTER
    grass.message(_("-- calcolo V --"))
    grass.run_command(
        "r.mapcalc", expression=f"num = {onecell}*v_i", overwrite=True
    )
    grass.run_command(
        "r.mapcalc", expression=f"den = {onecell}", overwrite=True
    )
    # TODO: parse correctly with python
    # numvar=`r.univar --q num|tail -n 1 |cut -f 2 -d " "`
    numvar = grass.parse_command("r.univar", map="num", quiet=True)
    # TODO: parse correctly with python
    # denvar=`r.univar --q den|tail -n 1 |cut -f 2 -d " "`
    denvar = grass.parse_command("r.univar", map="den", quiet=True)
    # TODO: check if digit after comma can processed without problems
    # Vfin=`echo "$numvar/$denvar"|bc -l`
    v_fin = numvar / denvar
    grass.message(_(f"Vfin: {v_fin}"))

    # TODO: check why there is an outcommented exit command here
    ####
    # exit
    ####

    #
    # calcolo I
    #
    grass.message(_("-- calcolo I -- ALT2 procedure"))

    grass.run_command(
        "v.db.addcolumn",
        map="su_segm",
        layer=2,
        columns=(
            "ai real, aj real, aall real, ci real, "
            "si real, cj real, sj real, num real"
        ),
    )
    grass.run_command(
        "db.execute",
        sql=(
            "update su_segm_edges_2 set ai = "
            "(select a_i from su_segm where cat=su_segm_edges_2.left)"
        ),
    )
    grass.run_command(
        "db.execute",
        sql=(
            "update su_segm_edges_2 set aj = "
            "(select a_i from su_segm where cat=su_segm_edges_2.right)"
        ),
    )
    grass.run_command(
        "db.execute",
        sql=(
            "update su_segm_edges_2 set aall = "
            "(select a_all from su_segm where cat=su_segm_edges_2.right)"
        ),
    )
    grass.run_command(
        "db.execute",
        sql=(
            "update su_segm_edges_2 set ci = "
            "cos(atan((sin(ai)+sin(aall))/(cos(ai)+cos(aall))))"
        ),
    )
    grass.run_command(
        "db.execute",
        sql=(
            "update su_segm_edges_2 set si = "
            "sin(atan((sin(ai)+sin(aall))/(cos(ai)+cos(aall))))"
        ),
    )
    grass.run_command(
        "db.execute",
        sql=(
            "update su_segm_edges_2 set cj = "
            "cos(atan((sin(aj)+sin(aall))/(cos(aj)+cos(aall))))"
        ),
    )
    grass.run_command(
        "db.execute",
        sql=(
            "update su_segm_edges_2 set sj = "
            "sin(atan((sin(aj)+sin(aall))/(cos(aj)+cos(aall))))"
        ),
    )
    grass.run_command(
        "db.execute",
        sql="update su_segm_edges_2 set num = ci*cj+si*sj",
    )

    # TODO: check if parsed correctly
    # num=`db.select -c sql=
    # "select sum(num) from su_segm_edges_2 where left<>-1 and right<>-1"`
    num = grass.parse_command(
        "db.select",
        flags="c",
        sql="select sum(num) from su_segm_edges_2 where left<>-1 and right<>-1",
    )
    # TODO: check if parsed correctly
    # den=`db.select -c sql=
    # "select count(*) from su_segm_edges_2 where left<>-1 and right<>-1"`
    den = grass.parse_command(
        "db.select",
        flags="c",
        sql="select count(*) from su_segm_edges_2 where left<>-1 and right<>-1",
    )
    grass.message(_(num))
    grass.message(_(den))

    # TODO: check if digit after comma can processed without problems
    # Ifin=`echo "$num/$den"|bc -l`
    i_fin = num / den
    grass.message(_(f"Ifin: {i_fin}"))

    # TODO: check why there is an outcommented exit command here
    ####
    # exit
    ####

    grass.run_command(
        "g.remove",
        type="raster",
        name=(
            "a_all,a_allrad,a_i,a_irad,areamap,aspect,base,basin_chk,cos,den,"
            "num,numero,sin,sumcos,sumcos_all,sumsin,sumsin_all,v_i,su_segm"
        ),
        flags="f",
    )
    grass.run_command(
        "g.remove", type="vector", name="su_segm,su_segm_edges", flags="f"
    )
    grass.run_command("r.mask", flags="r")

    # TODO: check if infos are appended correctly to file
    # echo "$areamin $cvmin $Vfin $Ifin" >> $outfile
    with open(outfile, "a") as file:
        file.write(f"{areamin} {cvmin} {v_fin} {i_fin}")
