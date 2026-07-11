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
import tempfile
import io
import os
from lammpsdoc import txt2html

class TestBasicFormatting(unittest.TestCase):
    def setUp(self):
        self.txt2html = txt2html.Txt2Html()

    def test_empty_string(self):
        self.assertEqual(self.txt2html.convert(""), "<HTML>\n"
                                                     "</HTML>\n")

    def test_single_paragraph(self):
        self.assertEqual(self.txt2html.convert("Hello World!\n"), "<HTML>\n"
                                                                   "<P>Hello World!\n"
                                                                   "</P>\n"
                                                                   "</HTML>\n")

    def test_two_paragraphs(self):
        s = self.txt2html.convert("Hello World!\n\nBye World!\n")
        self.assertEqual(s, "<HTML>\n"
                             "<P>Hello World!\n"
                             "</P>\n"
                             "<P>Bye World!\n"
                             "</P>\n"
                             "</HTML>\n")

    def test_line_concat(self):
        s = self.txt2html.convert("Hello World!\\\nBye World!\n")
        self.assertEqual(s, "<HTML>\n"
                             "<P>Hello World!Bye World!\n"
                             "</P>\n"
                             "</HTML>\n")

    def test_html_pass_through(self):
        s = self.txt2html.convert("<div>Raw HTML</div>\n")
        self.assertEqual(s, "<HTML>\n"
                             "<div>Raw HTML</div>\n\n"
                             "</HTML>\n")

    def test_ignore_rst(self):
        s = self.txt2html.convert("<!-- RST\n\n"
                                  ".. toctree::\n"
                                  "   :maxdepth: 2\n"
                                  "   :numbered:\n"
                                  "\n"
                                  "   Introduction\n"
                                  "END_RST -->\n")
        self.assertEqual("<HTML>\n"
                         "</HTML>\n", s)

    def test_ignore_html_only_markup(self):
        s = self.txt2html.convert("<!-- HTML_ONLY -->\n"
                                  "Hello World!\n"
                                  "<!-- END_HTML_ONLY -->\n")
        self.assertEqual("<HTML>\n"
                         "<!-- HTML_ONLY -->\n"
                         "Hello World!\n"
                         "<!-- END_HTML_ONLY -->\n\n"
                         "</HTML>\n", s)


class TestMarkup(unittest.TestCase):
    def setUp(self):
        self.markup = txt2html.HTMLMarkup()
        self.txt2html = txt2html.Txt2Html()

    def test_bold(self):
        self.assertEqual("<B>bold</B>", self.markup.convert("[bold]"))

    def test_italic(self):
        self.assertEqual("<I>italic</I>", self.markup.convert("{italic}"))

    def test_escape_markup(self):
        s = self.markup.convert("[bold] = \\[bold\\]\n"
                                "{italic} = \\{italic\\}\n")
        self.assertEqual("<B>bold</B> = [bold]\n"
                          "<I>italic</I> = {italic}\n", s)

    def test_link_markup(self):
        self.assertEqual("<A HREF = \"link\">Text</A>", self.markup.convert('"Text"_link'))

    def test_multiline_link_markup(self):
        s = self.txt2html.convert('"Te\n'
                                  'xt"_link')
        self.assertEqual("<HTML>\n"
                          "<P><A HREF = \"link\">Te\n"
                          "xt</A>\n"
                          "</P>\n"
                          "</HTML>\n", s)

    def test_ignore_punctuation_in_link(self):
        self.assertEqual("<A HREF = \"link\">Text</A>.", self.markup.convert('"Text"_link.'))
        self.assertEqual("<A HREF = \"link\">Text</A>,", self.markup.convert('"Text"_link,'))
        self.assertEqual("<A HREF = \"link\">Text</A>;", self.markup.convert('"Text"_link;'))
        self.assertEqual("<A HREF = \"link\">Text</A>:", self.markup.convert('"Text"_link:'))
        self.assertEqual("<A HREF = \"link\">Text</A>?", self.markup.convert('"Text"_link?'))
        self.assertEqual("<A HREF = \"link\">Text</A>!", self.markup.convert('"Text"_link!'))
        self.assertEqual("<A HREF = \"link\">Text</A>(", self.markup.convert('"Text"_link('))
        self.assertEqual("<A HREF = \"link\">Text</A>)", self.markup.convert('"Text"_link)'))

    def test_replace_alias_link(self):
        self.markup.add_link_alias("link", "replacement")
        self.assertEqual("<A HREF = \"replacement\">Text</A>", self.markup.convert('"Text"_link'))

class TestFormatting(unittest.TestCase):
    def setUp(self):
        self.txt2html = txt2html.Txt2Html()

    def test_paragraph_formatting(self):
        s = self.txt2html.convert("Hello :p\n")
        self.assertEqual(s, "<HTML>\n"
                             "<P>Hello \n"
                             "</P>\n"
                             "</HTML>\n")

    def test_two_paragraphs_through_formatting(self):
        text = "Hello :p\nBye :p\n"
        p = list(self.txt2html.paragraphs(text))
        s = self.txt2html.convert(text)
        self.assertEqual(len(p), 2)
        self.assertEqual(s, "<HTML>\n"
                             "<P>Hello \n"
                             "</P>\n"
                             "<P>Bye \n"
                             "</P>\n"
                             "</HTML>\n")

    def test_break_formatting(self):
        s = self.txt2html.convert("Hello :b\n")
        self.assertEqual(s, "<HTML>\n"
                             "Hello \n"
                             "<BR>\n"
                             "</HTML>\n")

    def test_preformat_formatting(self):
        s = self.txt2html.convert("Hello :pre\n")
        self.assertEqual(s, "<HTML>\n"
                             "<PRE>Hello \n"
                             "</PRE>\n"
                             "</HTML>\n")

    def test_center_formatting(self):
        s = self.txt2html.convert("Hello :c\n")
        self.assertEqual(s, "<HTML>\n"
                             "<CENTER>Hello \n"
                             "</CENTER>\n"
                             "</HTML>\n")

    def test_header_formatting(self):
        s = self.txt2html.convert("Level 1 :h1\n"
                                  "Level 2 :h2\n"
                                  "Level 3 :h3\n"
                                  "Level 4 :h4\n"
                                  "Level 5 :h5\n"
                                  "Level 6 :h6\n")
        self.assertEqual(s, "<HTML>\n"
                             "<H1>Level 1 \n"
                             "</H1>\n"
                             "<H2>Level 2 \n"
                             "</H2>\n"
                             "<H3>Level 3 \n"
                             "</H3>\n"
                             "<H4>Level 4 \n"
                             "</H4>\n"
                             "<H5>Level 5 \n"
                             "</H5>\n"
                             "<H6>Level 6 \n"
                             "</H6>\n"
                             "</HTML>\n")

    def test_all_paragraphs(self):
        s = self.txt2html.convert("one\n"
                                  "two\n"
                                  "three :all(p)\n")
        self.assertEqual("<HTML>\n"
                          "<P>one</P>\n"
                          "<P>two</P>\n"
                          "<P>three </P>\n"
                          "\n"
                          "</HTML>\n", s)

    def test_all_centered(self):
        s = self.txt2html.convert("one\n"
                                  "two\n"
                                  "three :all(c)\n")
        self.assertEqual("<HTML>\n"
                          "<CENTER>one</CENTER>\n"
                          "<CENTER>two</CENTER>\n"
                          "<CENTER>three </CENTER>\n"
                          "\n"
                          "</HTML>\n", s)

    def test_all_breaks(self):
        s = self.txt2html.convert("one\n"
                                  "two\n"
                                  "three :all(b)\n")
        self.assertEqual("<HTML>\n"
                          "one<BR>\n"
                          "two<BR>\n"
                          "three <BR>\n"
                          "\n"
                          "</HTML>\n", s)

    def test_links_with_all_breaks(self):
        s = self.txt2html.convert("\"one\"_link\n"
                                  "\"two\"_link\n"
                                  "\"three\"_link :all(b)\n")
        self.assertEqual("<HTML>\n"
                          "<A HREF = \"link\">one</A><BR>\n"
                          "<A HREF = \"link\">two</A><BR>\n"
                          "<A HREF = \"link\">three</A> <BR>\n"
                          "\n"
                          "</HTML>\n", s)

    def test_two_similar_links(self):
        s = self.txt2html.convert("\"one\"_linkA and \"one\"_linkAB\n")
        self.assertEqual("<HTML>\n"
                          "<P><A HREF = \"linkA\">one</A> and <A HREF = \"linkAB\">one</A>\n"
                          "</P>\n"
                          "</HTML>\n", s)

    def test_all_breaks_in_paragraph(self):
        s = self.txt2html.convert("one\n"
                                  "two\n"
                                  "three :all(b),p\n")
        self.assertEqual("<HTML>\n"
                          "<P>one<BR>\n"
                          "two<BR>\n"
                          "three <BR>\n"
                          "</P>\n"
                          "</HTML>\n", s)

    def test_all_list_items(self):
        s = self.txt2html.convert("one\n"
                                  "two\n"
                                  "three :all(l)\n")
        self.assertEqual("<HTML>\n"
                          "<LI>one\n"
                          "<LI>two\n"
                          "<LI>three \n"
                          "\n"
                          "</HTML>\n", s)

    def test_two_commands(self):
        s = self.txt2html.convert("one :ulb,l\n")
        self.assertEqual("<HTML>\n"
                          "<UL><LI>one \n"
                          "\n"
                          "</HTML>\n", s)

class TestListFormatting(unittest.TestCase):
    def setUp(self):
        self.txt2html = txt2html.Txt2Html()

    def test_unordered_list(self):
        s = self.txt2html.convert("one\n"
                                  "two\n"
                                  "three :ul\n")
        self.assertEqual(s, "<HTML>\n"
                             "<UL><LI>one\n"
                             "<LI>two\n"
                             "<LI>three \n"
                             "</UL>\n"
                             "</HTML>\n")

    def test_ordered_list(self):
        s = self.txt2html.convert("one\n"
                                  "two\n"
                                  "three :ol\n")
        self.assertEqual(s, "<HTML>\n"
                             "<OL><LI>one\n"
                             "<LI>two\n"
                             "<LI>three \n"
                             "</OL>\n"
                             "</HTML>\n")

    def test_elementwise_ordered_list(self):
        s = self.txt2html.convert("one :olb,l\n"
                                 "two :l\n"
                                 "three :ole,l\n")
        self.assertEqual("<HTML>\n"
                         "<OL><LI>one \n"
                         "\n"
                         "<LI>two \n"
                         "\n"
                         "<LI>three \n"
                         "</OL>\n"
                         "</HTML>\n", s)

    def test_definition_list(self):
        s = self.txt2html.convert("A\n"
                                  "first\n"
                                  "B\n"
                                  "second :dl\n")
        self.assertEqual(s, "<HTML>\n"
                             "<DL><DT>A\n"
                             "<DD>first\n"
                             "<DT>B\n"
                             "<DD>second \n"
                             "</DL>\n"
                             "</HTML>\n")

    def test_list_item(self):
        s = self.txt2html.convert("one :l\n")
        self.assertEqual(s, "<HTML>\n"
                             "<LI>one \n"
                             "\n"
                             "</HTML>\n")

    def test_definition_term(self):
        s = self.txt2html.convert("one :dt\n")
        self.assertEqual(s, "<HTML>\n"
                             "<DT>one \n"
                             "\n"
                             "</HTML>\n")

    def test_definition_description(self):
        s = self.txt2html.convert("one :dd\n")
        self.assertEqual(s, "<HTML>\n"
                             "<DD>one \n"
                             "\n"
                             "</HTML>\n")

    def test_unordered_list_begin(self):
        s = self.txt2html.convert("one :ulb\n")
        self.assertEqual(s, "<HTML>\n"
                             "<UL>one \n"
                             "\n"
                             "</HTML>\n")

    def test_unordered_list_end(self):
        s = self.txt2html.convert("one :ule\n")
        self.assertEqual(s, "<HTML>\n"
                             "one \n"
                             "</UL>\n"
                             "</HTML>\n")

    def test_ordered_list_begin(self):
        s = self.txt2html.convert("one :olb\n")
        self.assertEqual(s, "<HTML>\n"
                             "<OL>one \n"
                             "\n"
                             "</HTML>\n")

    def test_ordered_list_end(self):
        s = self.txt2html.convert("one :ole\n")
        self.assertEqual(s, "<HTML>\n"
                             "one \n"
                             "</OL>\n"
                             "</HTML>\n")

    def test_definition_list_begin(self):
        s = self.txt2html.convert("one :dlb\n")
        self.assertEqual(s, "<HTML>\n"
                             "<DL>one \n"
                             "\n"
                             "</HTML>\n")

    def test_definition_list_end(self):
        s = self.txt2html.convert("one :dle\n")
        self.assertEqual(s, "<HTML>\n"
                             "one \n"
                             "</DL>\n"
                             "</HTML>\n")

class TestSpecialCommands(unittest.TestCase):
    def setUp(self):
        self.txt2html = txt2html.Txt2Html()

    def test_line(self):
        s = self.txt2html.convert("one :line\n")
        self.assertEqual(s, "<HTML>\n"
                             "<HR>one \n"
                             "\n"
                             "</HTML>\n")

    def test_image(self):
        s = self.txt2html.convert("one :image(file)\n")
        self.assertEqual(s, "<HTML>\n"
                             "<IMG SRC = \"file\">one \n"
                             "\n"
                             "</HTML>\n")

    def test_image_with_link(self):
        s = self.txt2html.convert("one :image(file,link)\n")
        self.assertEqual(s, "<HTML>\n"
                             "<A HREF = \"link\"><IMG SRC = \"file\"></A>one \n"
                             "\n"
                             "</HTML>\n")

    def test_named_link(self):
        s = self.txt2html.convert("one :link(name)\n")
        self.assertEqual(s, "<HTML>\n"
                             "<A NAME = \"name\"></A>one \n"
                             "\n"
                             "</HTML>\n")

    def test_define_link_alias(self):
        s = self.txt2html.convert("one :link(alias,value)\n"
                                  "\"test\"_alias")
        self.assertEqual(s, "<HTML>\n"
                             "one \n"
                             "\n"
                             "<P><A HREF = \"value\">test</A>\n"
                             "</P>\n"
                             "</HTML>\n")

    def test_define_link_alias_later(self):
        s = self.txt2html.convert("\"test\"_alias\n\n"
                                  "one :link(alias,value)\n")
        self.assertEqual(s, "<HTML>\n"
                             "<P><A HREF = \"value\">test</A>\n"
                             "</P>\n"
                             "one \n"
                             "\n"
                             "</HTML>\n")

class TestTableCommand(unittest.TestCase):
    def setUp(self):
        self.txt2html = txt2html.Txt2Html()

    def test_single_table_row(self):
        s = self.txt2html.convert("a,b,c :tb")
        self.assertEqual("<HTML>\n"
                          "<DIV ALIGN=center><TABLE  BORDER=1 >\n"
                          "<TR><TD >a</TD><TD >b</TD><TD >c \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "</HTML>\n", s)

    def test_multiple_table_rows(self):
        s = self.txt2html.convert("a,b,c,\nd,e,f,\ng,h,i :tb")
        self.assertEqual("<HTML>\n"
                          "<DIV ALIGN=center><TABLE  BORDER=1 >\n"
                          "<TR><TD >a</TD><TD >b</TD><TD >c</TD><TD ></TD></TR>\n"
                          "<TR><TD >d</TD><TD >e</TD><TD >f</TD><TD ></TD></TR>\n"
                          "<TR><TD >g</TD><TD >h</TD><TD >i \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "</HTML>\n", s)

    def test_fixed_table_columns(self):
        s = self.txt2html.convert("a,b,c,d,\ne,f,g,h :tb(c=3)\n")
        self.assertEqual("<HTML>\n"
                          "<DIV ALIGN=center><TABLE  BORDER=1 >\n"
                          "<TR><TD >a</TD><TD >b</TD><TD >c</TD></TR>\n"
                          "<TR><TD >d</TD><TD >e</TD><TD >f</TD></TR>\n"
                          "<TR><TD >g</TD><TD >h \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "</HTML>\n", s)

    def test_change_cell_separator(self):
        s = self.txt2html.convert("a:b:c :tb(s=:)\n")
        self.assertEqual("<HTML>\n"
                          "<DIV ALIGN=center><TABLE  BORDER=1 >\n"
                          "<TR><TD >a</TD><TD >b</TD><TD >c \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "</HTML>\n", s)

    def test_change_table_border(self):
        s = self.txt2html.convert("a,b,c :tb(b=0)")
        self.assertEqual("<HTML>\n"
                          "<DIV ALIGN=center><TABLE  BORDER=0 >\n"
                          "<TR><TD >a</TD><TD >b</TD><TD >c \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "</HTML>\n", s)

    def test_change_cell_width(self):
        s = self.txt2html.convert("a,b,c :tb(w=10)")
        self.assertEqual("<HTML>\n"
                          "<DIV ALIGN=center><TABLE  BORDER=1 >\n"
                          "<TR><TD WIDTH=\"10\">a</TD><TD WIDTH=\"10\">b</TD><TD WIDTH=\"10\">c \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "</HTML>\n", s)

    def test_change_table_width(self):
        s = self.txt2html.convert("a,b,c :tb(w=10%)")
        self.assertEqual("<HTML>\n"
                          "<DIV ALIGN=center><TABLE  WIDTH=\"10%\" BORDER=1 >\n"
                          "<TR><TD >a</TD><TD >b</TD><TD >c \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "</HTML>\n", s)

    def test_change_table_alignment(self):
        s = self.txt2html.convert("a :tb(a=l)\n"
                                  "b :tb(a=c)\n"
                                  "c :tb(a=r)\n")
        self.assertEqual("<HTML>\n"
                          "<DIV ALIGN=left><TABLE  BORDER=1 >\n"
                          "<TR><TD >a \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "<DIV ALIGN=center><TABLE  BORDER=1 >\n"
                          "<TR><TD >b \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "<DIV ALIGN=right><TABLE  BORDER=1 >\n"
                          "<TR><TD >c \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "</HTML>\n", s)

    def test_change_cell_alignment(self):
        s = self.txt2html.convert("a :tb(ea=l)\n"
                                  "b :tb(ea=c)\n"
                                  "c :tb(ea=r)\n")
        self.assertEqual("<HTML>\n"
                          "<DIV ALIGN=center><TABLE  BORDER=1 >\n"
                          "<TR ALIGN=\"left\"><TD >a \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "<DIV ALIGN=center><TABLE  BORDER=1 >\n"
                          "<TR ALIGN=\"center\"><TD >b \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "<DIV ALIGN=center><TABLE  BORDER=1 >\n"
                          "<TR ALIGN=\"right\"><TD >c \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "</HTML>\n", s)

    def test_change_cell_vertical_alignment(self):
        s = self.txt2html.convert("a :tb(eva=t,ea=l)\n"
                                  "b :tb(eva=m)\n"
                                  "c :tb(eva=ba)\n"
                                  "d :tb(eva=bo)\n")
        self.assertEqual("<HTML>\n"
                          "<DIV ALIGN=center><TABLE  BORDER=1 >\n"
                          "<TR ALIGN=\"left\" VALIGN =\"top\"><TD >a \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "<DIV ALIGN=center><TABLE  BORDER=1 >\n"
                          "<TR VALIGN =\"middle\"><TD >b \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "<DIV ALIGN=center><TABLE  BORDER=1 >\n"
                          "<TR VALIGN =\"baseline\"><TD >c \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "<DIV ALIGN=center><TABLE  BORDER=1 >\n"
                          "<TR VALIGN =\"bottom\"><TD >d \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "</HTML>\n", s)

    def test_custom_column_width(self):
        s = self.txt2html.convert("a,b,c :tb(cw1=30)\n"
                                  "a,b,c :tb(cw2=30%)\n")
        self.assertEqual("<HTML>\n"
                          "<DIV ALIGN=center><TABLE  BORDER=1 >\n"
                          "<TR><TD WIDTH=\"30\">a</TD><TD >b</TD><TD >c \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "<DIV ALIGN=center><TABLE  BORDER=1 >\n"
                          "<TR><TD >a</TD><TD WIDTH=\"30%\">b</TD><TD >c \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "</HTML>\n", s)

    def test_custom_column_alignment(self):
        s = self.txt2html.convert("a,b,c :tb(w=30,ca1=l)\n"
                                  "a,b,c :tb(ca2=c)\n"
                                  "a,b,c :tb(ca3=r)\n")
        self.assertEqual("<HTML>\n"
                          "<DIV ALIGN=center><TABLE  BORDER=1 >\n"
                          "<TR><TD WIDTH=\"30\" ALIGN =\"left\">a</TD><TD WIDTH=\"30\">b</TD><TD WIDTH=\"30\">c \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "<DIV ALIGN=center><TABLE  BORDER=1 >\n"
                          "<TR><TD >a</TD><TD  ALIGN =\"center\">b</TD><TD >c \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "<DIV ALIGN=center><TABLE  BORDER=1 >\n"
                          "<TR><TD >a</TD><TD >b</TD><TD  ALIGN =\"right\">c \n"
                          "</TD></TR></TABLE></DIV>\n\n"
                          "</HTML>\n", s)

class TestTxt2HtmlCLI(unittest.TestCase):
    def setUp(self):
        self.out = io.StringIO()
        self.err = io.StringIO()
        self.app = txt2html.Txt2HtmlConverter()

    def test_convert_single_file(self):
        with tempfile.NamedTemporaryFile(mode='w+t') as f:
            f.write('Hello World!\n')
            f.flush()
            args = [f.name]
            self.app.run(args=args, out=self.out, err=self.err)
            self.assertEqual("<HTML>\n"
                              "<P>Hello World!\n"
                              "</P>\n"
                              "</HTML>\n", self.out.getvalue())
            self.assertEqual("Converting " + f.name + " ...\n", self.err.getvalue())

    def test_convert_multiple_files(self):
        with tempfile.NamedTemporaryFile(mode='w+t') as f:
            with tempfile.NamedTemporaryFile(mode='w+t') as g:
                f.write('Hello World!\n')
                f.flush()
                g.write('Hello World!\n')
                g.flush()
                args = [f.name, g.name]
                self.app.run(args=args, out=self.out, err=self.err)
                self.assertEqual("", self.out.getvalue())
                self.assertEqual("Converting " + f.name + " ...\n"
                                  "Converting " + g.name + " ...\n", self.err.getvalue())
                self.assertTrue(os.path.exists(f.name + ".html"))
                self.assertTrue(os.path.exists(g.name + ".html"))
                os.remove(f.name + ".html")
                os.remove(g.name + ".html")

    def test_break_flag(self):
        with tempfile.NamedTemporaryFile(mode='w+t') as f:
            f.write('Hello World!\n')
            f.flush()
            args = ["-b", f.name]
            self.app.run(args=args, out=self.out, err=self.err)
            self.assertEqual("<HTML>\n"
                              "<P>Hello World!\n"
                              "</P>\n"
                              "<!-- PAGE BREAK -->\n"
                              "</HTML>\n", self.out.getvalue())

    def test_skip_files(self):
        with tempfile.NamedTemporaryFile(mode='w+t') as f:
            with tempfile.NamedTemporaryFile(mode='w+t') as g:
                f.write('Hello World!\n')
                f.flush()
                g.write('Hello World!\n')
                g.flush()
                args = ["-x", g.name, f.name, g.name]
                self.app.run(args=args, out=self.out, err=self.err)
                self.assertEqual("", self.out.getvalue())
                self.assertEqual("Converting " + f.name + " ...\n", self.err.getvalue())
                self.assertTrue(os.path.exists(f.name + ".html"))
                self.assertFalse(os.path.exists(g.name + ".html"))
                os.remove(f.name + ".html")

    def test_create_title_option(self):
        with tempfile.NamedTemporaryFile(mode='w+t') as f:
            f.write('Hello World! :h1\n')
            f.flush()
            args = ["--generate-title", f.name]
            self.app.run(args=args, out=self.out, err=self.err)
            self.assertEqual("<HTML>\n"
                              "<HEAD>\n"
                              "<TITLE>Hello World!</TITLE>\n"
                              "</HEAD>\n"
                              "<H1>Hello World! \n"
                              "</H1>\n"
                              "</HTML>\n", self.out.getvalue())

class TestMathMarkup(unittest.TestCase):
    def setUp(self):
        self.txt2html = txt2html.Txt2Html()

    def test_detect_latex_equation(self):
        s = self.txt2html.convert("\\begin\\{equation\\} T_\\{ij\\}(r_\\{ij\\}) = 1 - \\left( 1 +\n"
                                  "\\frac\{s_\\{ij\} r_\\{ij\} \\}\\{2\\} \\right)\n"
                                  "\\exp \\left( - s_\\{ij\\} r_\\{ij\\} \\right) \\end\\{equation\\}\n")
        self.assertEqual("<HTML>\n"
                         "<P>\\begin{equation} T_{ij}(r_{ij}) = 1 - \\left( 1 +\n"
                         "\\frac{s_{ij} r_{ij} }{2} \\right)\n"
                         "\\exp \\left( - s_{ij} r_{ij} \\right) \\end{equation}\n"
                         "</P>\n"
                         "</HTML>\n", s)

if __name__ == '__main__':
    unittest.main()
