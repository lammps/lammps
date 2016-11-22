"""Tests that the Lennard-Jones test case works properly.

Copyright (C) 2013, Joshua More and Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.
"""

from common import TestSimulation

def test_lj_gas():
    ts = TestSimulation(input="../../test/lj/gas/input.xml", driver="../../drivers/driver.x")
    ts.run()
    # Test properties (e.g. latest positions/temperature etc)
