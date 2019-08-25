# LAMMPS Documentation Utilities
#
# Copyright (C) 2015 Richard Berger
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import unittest
from lammpsdoc import txt2rst

class TestStructuralFilters(unittest.TestCase):
    def setUp(self):
        self.txt2rst = txt2rst.Txt2Rst()

    def test_detect_and_replace_warnings(self):
        s = self.txt2rst.convert("IMPORTANT NOTE: Content\n")
        self.assertEqual(".. warning::\n\n"
                         "   Content\n"
                         "\n", s)

    def test_detect_and_replace_note(self):
        s = self.txt2rst.convert("NOTE: Content\n")
        self.assertEqual(".. note::\n\n"
                         "   Content\n"
                         "\n", s)

    def test_detect_command_and_add_to_index(self):
        s = self.txt2rst.convert("some command\n")
        self.assertEqual(".. index:: some\n\n"
                         "some command\n\n", s)

    def test_filter_file_header_and_append_common_links(self):
        s = self.txt2rst.convert("some random text\n"
                                 "which should be ignored\n"
                                 "----------\n\n"
                                 "Title\n")
        self.assertEqual("Title\n\n"
                         "\n.. _lws: http://lammps.sandia.gov\n"
                         ".. _ld: Manual.html\n"
                         ".. _lc: Commands_all.html\n", s)

    def test_filter_multiple_horizontal_rules(self):
        s = self.txt2rst.convert(":hline\n"
                                 "   \n\n"
                                 ":hline\n")
        self.assertEqual("\n\n", s)
