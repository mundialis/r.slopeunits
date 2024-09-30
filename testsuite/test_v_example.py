#!/usr/bin/env python3
"""
############################################################################
#
# MODULE:      v.tiles test
# AUTHOR(S):   Lina Krisztian

# PURPOSE:      Tests v.tiles
# COPYRIGHT:   (C) 2022 by mundialis GmbH & Co. KG and the GRASS Development
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

import os

from grass.gunittest.case import TestCase
from grass.gunittest.main import test
from grass.gunittest.gmodules import SimpleModule
import grass.script as grass


class TestVTiles(TestCase):
    """Test class with as an example for the example addon v.example."""

    pid = os.getpid()
    rm_vec_class = []
    rm_vec_test = []
    region_vec_file = os.path.join("data", "area_beuel.geojson")
    region_vec = f"beuel_vec_{pid}"
    rm_vec_class.append(region_vec)
    output_base = f"tile_{pid}"
    rm_vec_test.append(output_base)

    @classmethod
    # pylint: disable=invalid-name
    def setUpClass(cls):
        """Ensures expected computational region and generated data"""
        # import and set region
        cls.runModule(
            "v.import", input=cls.region_vec_file, output=cls.region_vec
        )
        cls.runModule("g.region", vector=cls.region_vec, res=1000, flags="a")

    @classmethod
    # pylint: disable=invalid-name
    def tearDownClass(cls):
        """Remove the temporary region and generated data"""
        for vec in cls.rm_vec_class:
            cls.runModule("g.remove", type="vector", name=vec, flags="f")

    # pylint: disable=invalid-name
    def tearDown(self):
        """Remove the outputs created
        This is executed after each test run.
        """
        for vec in self.rm_vec_test:
            self.runModule(
                "g.remove", type="vector", pattern=f"{vec}*", flags="f"
            )

    def test_v_tiles(self):
        """Test v.tiles without polygon_aoi map given"""
        v_check = SimpleModule(
            "v.tiles", output=self.output_base, box=[1000, 1000]
        )
        self.assertModule(v_check)
        # check if output is the desired one
        tiles = grass.parse_command(
            "g.list", pattern=f"{self.output_base}*", type="vector"
        )
        # expect creation of 56 tiles
        self.assertTrue(len(tiles) == 56)

    def test_v_tiles_with_polygonaoi(self):
        """Test v.tiles with polygon_aoi map given"""
        v_check = SimpleModule(
            "v.tiles",
            output=self.output_base,
            box=[1000, 1000],
            polygon_aoi=self.region_vec,
        )
        self.assertModule(v_check)
        # check if output is the desired one
        tiles = grass.parse_command(
            "g.list", pattern=f"{self.output_base}*", type="vector"
        )
        # expect creation of 37 tiles
        self.assertTrue(len(tiles) == 37)


if __name__ == "__main__":
    test()
