#!/usr/bin/env python3
"""
############################################################################
#
# MODULE:      v.example
# AUTHOR(S):   {NAME}

# PURPOSE:     {SHORT DESCRIPTION}
# COPYRIGHT:   (C) {YEAR} by mundialis GmbH & Co. KG and the GRASS Development
#              Team
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
#############################################################################
"""
# %Module
# % description: {SHORT DESCRIPTION}.
# % keyword: vector
# % keyword: grid
# %end

# %option G_OPT_V_OUTPUT
# % key: output
# % required: yes
# % description: Base name for tiles.
# %end

# %option
# % key: box
# % key_desc: width,height
# % type: double
# % required: yes
# % multiple: no
# % description: Width and height of boxes in grid
# %end

# %option G_OPT_V_INPUT
# % key: polygon_aoi
# % required: no
# % description: Polygon to extract grid-tiles for
# %end

# %option G_OPT_V_FIELD
# %end

# %option G_OPT_DB_COLUMN
# %end

# %flag
# % key: p
# % description: Print attribute values of given column for input vector polygon_aoi
# %end


# import needed libraries
import atexit
import os
import grass.script as grass

# initialize global variables
rm_vec = []


# cleanup function (can be extended)
def cleanup():
    """Cleanup fuction (can be extended)"""
    nulldev = open(os.devnull, "w", encoding="utf-8")
    kwargs = {"flags": "f", "quiet": True, "stderr": nulldev}
    for rmvec in rm_vec:
        if grass.find_file(name=rmvec, element="vector")["file"]:
            grass.run_command("g.remove", type="vector", name=rmvec, **kwargs)


def main():
    """Main function of v.example"""
    global rm_vec

    # print attribute values if requested
    if flags["p"] and options["polygon_aoi"] and options["column"]:
        aoi_vector = options["polygon_aoi"]
        layer = options["layer"]
        column = options["column"]

        # examples for interfaces in
        # python/grass/script/vector.py
        # python/grass/script/db.py

        # get all database connections as dictionary
        db_connections = grass.vector_db(aoi_vector)
        if len(db_connections) == 0:
            grass.fatal(_("No database connections for map %s") % aoi_vector)

        # check if there is an attribute table at the given layer
        # grass.vector_layer_db() will exit with fatal error if there is
        # no database connection for the given layer
        db_connection = grass.vector_layer_db(aoi_vector, layer)

        # get attribute columns as dictionary (default)
        columns = grass.vector_columns(aoi_vector, layer, getDict=True)
        if column not in columns.keys():
            grass.fatal(
                _("Column %s does not exist in layer %s of vector %s")
                % (column, layer, aoi_vector)
            )

        # get attribute columns as list
        columns = grass.vector_columns(aoi_vector, layer, getDict=False)
        colidx = -1
        try:
            colidx = columns.index(column)
        except ValueError:
            grass.fatal(
                _("Column %s does not exist in layer %s of vector %s")
                % (column, layer, aoi_vector)
            )

        if colidx >= 0:
            grass.verbose(
                _("Found column %s in vector %s, layer %s")
                % (column, aoi_vector, layer)
            )

        # query the database and table directly using information in the
        # database connection of the vector
        database = db_connection["database"]
        table = db_connection["table"]
        driver = db_connection["driver"]

        table_description = grass.db_describe(
            table=table, database=database, driver=driver
        )
        found = False
        # TODO: pythonize this for loop
        for i in range(len(table_description["cols"])):
            column_description = table_description["cols"][i]
            if column_description[0] == column:
                found = True
                break

        if found is False:
            grass.fatal(
                _("Column %s does not exist in layer %s of vector %s")
                % (column, layer, aoi_vector)
            )

        # select attribute values with vector_db_select()
        grass.message(
            _("Print attribute values using %s") % "vector_db_select()"
        )
        column_values = grass.vector_db_select(
            aoi_vector, int(layer), columns=column
        )
        # go over table rows
        for key in column_values["values"]:
            # print value of first selected column
            print(column_values["values"][key][0])

        # select attribute values with SQL statement
        grass.message(_("Print attribute values using %s") % "db_select()")
        sql = f"select {column} from {table}"
        values = grass.db_select(sql=sql, database=database, driver=driver)
        for value in values:
            print(value[0])

    # use process ID as suffix for unique names of temporary maps
    pid = os.getpid()

    # set grid
    out_grid = f"kacheln_{pid}"
    rm_vec.append(out_grid)
    grass.run_command(
        "v.mkgrid",
        map=out_grid,
        position="region",
        box=options["box"],
        quiet=True,
    )
    # extract only polygon_aoi area if given:
    if options["polygon_aoi"]:
        out_overlay = f"overlay_aoi_grid_{pid}"
        rm_vec.append(out_overlay)
        grass.run_command(
            "v.overlay",
            ainput=out_grid,
            binput=options["polygon_aoi"],
            operator="and",
            output=out_overlay,
        )
    else:
        out_overlay = out_grid
    # divide into tiles
    kachel_num = grass.parse_command(
        "v.db.select", map=out_overlay, columns="cat", flags="c", quiet=True
    )
    for kachel in kachel_num:
        grass.run_command(
            "v.extract",
            input=out_overlay,
            output=f"{options['output']}_{kachel}",
            cats=kachel,
            quiet=True,
        )
    grass.message(_(f"Created {len(kachel_num)} tiles."))


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
