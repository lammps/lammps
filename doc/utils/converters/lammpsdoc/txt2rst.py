#! /usr/bin/env python3
# LAMMPS Documentation Utilities
#
# Converter of LAMMPS documentation format to Sphinx ReStructured Text
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

import os
import re
import argparse
from lammpsdoc import lammps_filters
from lammpsdoc.txt2html import Markup, Formatting, TxtParser, TxtConverter


class RSTMarkup(Markup):
    def __init__(self):
        super().__init__()

    def bold_start(self):
        return "**"

    def bold_end(self):
        return "**"

    def italic_start(self):
        return "*"

    def italic_end(self):
        return "*"

    def bold(self, text):
        """ RST requires a space after inline formats.
        For words which only partially apply a format add a backslash and whitespace to create valid RST"""
        text = re.sub(r'([^\s\\])\[([^\]\\]+)\]', r'\1\\ [\2]', text)
        text = re.sub(r'([^\\]?)\[([^\]\\]+)\]([^\s])', r'\1[\2]\\ \3', text)
        text = super().bold(text)
        return text

    def italic(self, text):
        """ RST requires a space after inline formats.
        For words which only partially apply a format add a backslash and whitespace to create valid RST"""
        text = re.sub(r'([^\s\\])\{([^\}\\]+)\}', r'\1\\ {\2}', text)
        text = re.sub(r'([^\\]?)\{([^\}\\]+)\}([^\s])', r'\1{\2}\\ \3', text)
        text = super().italic(text)
        return text

    def convert(self, text):
        text = self.escape_rst_chars(text)
        text = super().convert(text)
        text = self.inline_math(text)
        return text

    def escape_rst_chars(self, text):
        text = text.replace('*', '\\*')
        text = text.replace('^', '\\^')
        text = text.replace('|', '\\|')
        text = re.sub(r'([^"])_([ \t\n\r\f])', r'\1\\\\_\2', text)
        text = re.sub(r'([^"])_([^ \t\n\r\f])', r'\1\\_\2', text)
        return text

    def unescape_rst_chars(self, text):
        text = text.replace('\\*', '*')
        text = text.replace('\\^', '^')
        text = self.unescape_underscore(text)
        text = text.replace('\\|', '|')
        return text

    def unescape_underscore(self, text):
        return text.replace('\\_', '_')

    def inline_math(self, text):
        start_pos = text.find("\\(")
        end_pos = text.find("\\)")

        while start_pos >= 0 and end_pos >= 0:
            original = text[start_pos:end_pos+2]
            formula = original[2:-2]
            formula = self.unescape_rst_chars(formula)
            replacement = ":math:`" + formula.replace('\n', ' ').strip() + "`"
            text = text.replace(original, replacement)

            start_pos = text.find("\\(")
            end_pos = text.find("\\)")

        return text

    def create_link(self, content, href):
        content = content.strip()
        content = content.replace('\n', ' ')

        href = self.unescape_rst_chars(href)

        anchor_pos = href.find('#')

        if anchor_pos >= 0 and not href.startswith('http'):
            href = href[anchor_pos+1:]
            return ":ref:`%s <%s>`" % (content, href)

        if href in self.references:
            return ":ref:`%s <%s>`" % (content, href)
        elif href in self.aliases:
            href = "%s_" % href
        elif href.endswith('.html') and not href.startswith('http') and 'USER/atc' not in href:
            href = href[0:-5]
            return ":doc:`%s <%s>`" % (content, href)

        return "`%s <%s>`_" % (content, href)


class RSTFormatting(Formatting):
    RST_HEADER_TYPES = '#*=-^"'

    def __init__(self, markup):
        super().__init__(markup)
        self.indent_level = 0

    def paragraph(self, content):
        if self.indent_level > 0:
            return '\n' + self.list_indent(content.strip(), self.indent_level)

        return content.strip() + "\n"

    def center(self, content):
        return content

    def linebreak(self, content):
        return content.strip()

    def preformat(self, content):
        content = self.markup.unescape_underscore(content)
        if self.indent_level > 0:
            return self.list_indent("\n.. parsed-literal::\n\n" + self.indent(content.rstrip()), self.indent_level)
        return "\n.. parsed-literal::\n\n" + self.indent(content.rstrip())

    def horizontal_rule(self, content):
        return "\n----------\n\n" + content.strip()

    def image(self, content, file, link=None):
        # 2017-12-07: commented out to disable thumbnail processing due to dropping
        #             support for obsolete sphinxcontrib.images extension
        #
        #if link and (link.lower().endswith('.jpg') or
        #                 link.lower().endswith('.jpeg') or
        #                 link.lower().endswith('.png') or
        #                 link.lower().endswith('.gif')):
        #    converted = ".. thumbnail:: " + self.markup.unescape_rst_chars(link) + "\n"
        #else:
        converted = ".. image:: " + self.markup.unescape_rst_chars(file) + "\n"
        if link:
            converted += "   :target: " + self.markup.unescape_rst_chars(link) + "\n"

        if "c" in self.current_command_list:
            converted += "   :align: center\n"

        return converted + content.strip()

    def named_link(self, paragraph, name):
        self.markup.add_internal_reference(name)
        return (".. _%s:\n\n" % name) + paragraph

    def define_link_alias(self, paragraph, alias, value):
        self.markup.add_link_alias(alias, value)
        return (".. _%s: %s\n\n" % (alias, value)) + paragraph

    def header(self, content, level):
        header_content = content.strip()
        header_content = re.sub(r'[0-9]+\.([0-9]*\.?)*\s+', '', header_content)
        header_underline = RSTFormatting.RST_HEADER_TYPES[level-1] * len(header_content)
        return header_content + "\n" + header_underline + "\n"

    def unordered_list_item(self, paragraph):
        return "* " + paragraph.strip().replace('\n', '\n  ')

    def ordered_list_item(self, paragraph, index):
        if index is None:
            index = "#"
        return str(index) + ". " + paragraph.strip().replace('\n', '\n   ')

    def definition_term(self, paragraph):
        return paragraph.strip()

    def definition_description(self, paragraph):
        return self.indent(paragraph.strip())

    def unordered_list_begin(self, paragraph):
        self.indent_level += 1
        return paragraph

    def unordered_list_end(self, paragraph):
        self.indent_level -= 1
        return paragraph.rstrip() + '\n'

    def ordered_list_begin(self, paragraph):
        if paragraph.startswith('* '):
            paragraph = '#. ' + paragraph[2:]
        self.indent_level += 1
        return paragraph

    def definition_list_begin(self, paragraph):
        return paragraph

    def definition_list_end(self, paragraph):
        return paragraph

    def ordered_list_end(self, paragraph):
        self.indent_level -= 1
        return paragraph.rstrip() + '\n'

    def ordered_list(self, paragraph):
        paragraph = super().ordered_list(paragraph)
        return paragraph.rstrip() + '\n'

    def unordered_list(self, paragraph):
        paragraph = super().unordered_list(paragraph)
        return paragraph.rstrip() + '\n'

    def all_breaks(self, paragraph):
        indented = ""
        for line in paragraph.splitlines():
            indented += "| %s\n" % line
        indented += "| \n"
        return indented

    def begin_document(self):
        return ""

    def end_document(self):
        return ""

    def raw_html(self, content):
        raw_directive = ".. raw:: html\n\n"
        return raw_directive + self.indent(content)

    def indent(self, content):
        indented = ""
        for line in content.splitlines():
            indented += "   %s\n" % line
        return indented

    def list_indent(self, content, level=1):
        indented = ""
        for line in content.splitlines():
            indented += "  " * level + ("%s\n" % line)
        return indented

    def get_max_column_widths(self, rows):
        num_columns = max([len(row) for row in rows])
        max_widths = [0] * num_columns

        for columns in rows:
            for col_idx, column in enumerate(columns):
                max_widths[col_idx] = max(max_widths[col_idx], len(column.strip())+2)

        return max_widths

    def create_table_horizontal_line(self, max_widths):
        cell_borders = ['-' * width for width in max_widths]
        return '+' + '+'.join(cell_borders) + "+"

    def table(self, paragraph, configuration):
        paragraph = self.protect_rst_directives(paragraph)

        if configuration['num_columns'] == 0:
            rows = self.create_table_with_columns_based_on_newlines(paragraph, configuration['separator'])
        else:
            rows = self.create_table_with_fixed_number_of_columns(paragraph, configuration['separator'],
                                                                  configuration['num_columns'])

        column_widths = self.get_max_column_widths(rows)
        max_columns = len(column_widths)
        horizontal_line = self.create_table_horizontal_line(column_widths) + "\n"

        tbl = horizontal_line

        for row_idx in range(len(rows)):
            columns = rows[row_idx]

            tbl += "| "
            for col_idx in range(max_columns):
                if col_idx < len(columns):
                    col = columns[col_idx].strip()
                else:
                    col = ""

                tbl += col.ljust(column_widths[col_idx]-2, ' ')
                tbl += " |"

                if col_idx < max_columns - 1:
                    tbl += " "
            tbl += "\n"
            tbl += horizontal_line

        tbl = self.restore_rst_directives(tbl)
        return tbl

    def protect_rst_directives(self, content):
        content = content.replace(":doc:", "0DOC0")
        content = content.replace(":ref:", "0REF0")
        return content

    def restore_rst_directives(self, content):
        content = content.replace("0DOC0", ":doc:")
        content = content.replace("0REF0", ":ref:")
        return content

    def math(self, content):
        eqs = content.split(r'\end{equation}')

        text = ""

        if len(eqs) > 1:
            post = eqs[-1].strip()
            eqs = eqs[0:-1]
        else:
            post = ""

        for eq in eqs:
            if len(eq.strip()) == 0:
                continue
            parts = eq.split(r'\begin{equation}')
            assert(len(parts) == 1 or len(parts) == 2)
            if len(parts) == 2:
                start = parts[0].strip()
                body = parts[1]
            else:
                start = ""
                body = parts[0]

            body = self.markup.unescape_rst_chars(body)

            if len(start) > 0:
                text += start + "\n"
            text += "\n.. math::\n\n"
            text += self.indent(r'\begin{equation}' + body.strip('\n') + r'\end{equation}')
            text += "\n"

        return text + post


class Txt2Rst(TxtParser):
    def __init__(self):
        super().__init__()
        self.markup = RSTMarkup()
        self.format = RSTFormatting(self.markup)
        self.register_filters()

    def register_filters(self):
        self.paragraph_filters.append(lammps_filters.detect_and_format_notes)
        self.document_filters.append(lammps_filters.filter_file_header_until_first_horizontal_line)
        self.document_filters.append(lammps_filters.detect_and_add_command_to_index)
        self.document_filters.append(lammps_filters.filter_multiple_horizontal_rules)
        self.document_filters.append(lammps_filters.promote_doc_keywords)
        self.document_filters.append(lammps_filters.merge_preformatted_sections)

    def is_ignored_textblock_begin(self, line):
        return line.startswith('<!-- HTML_ONLY -->')

    def is_ignored_textblock_end(self, line):
        return line.startswith('<!-- END_HTML_ONLY -->')

    def is_raw_textblock_begin(self, line):
        return line.startswith('<!-- RST')

    def is_raw_textblock_end(self, line):
        return line.startswith('END_RST -->')

    def order_commands(self, commands):
        if 'ule' in commands and 'l' in commands and commands.index('ule') >  commands.index('l'):
            return commands
        elif 'ole' in commands and 'l' in commands and commands.index('ole') > commands.index('l'):
            return commands
        return super().order_commands(commands)

    def transform_paragraphs(self, content):
        if self.format.indent_level > 0:
            raise Exception("unbalanced number of ulb,ule or olb,ole pairs!")
        return super().transform_paragraphs(content)


class Txt2RstConverter(TxtConverter):
    def get_argument_parser(self):
        parser = argparse.ArgumentParser(description='converts a text file with simple formatting & markup into '
                                                     'Restructured Text for Sphinx.')
        parser.add_argument('-x', metavar='file-to-skip', dest='skip_files', action='append')
        parser.add_argument('--verbose', '-v', dest='verbose', action='store_true')
        parser.add_argument('--output-directory', '-o', dest='output_dir')
        parser.add_argument('files',  metavar='file', nargs='+', help='one or more files to convert')
        return parser

    def create_converter(self, args):
        return Txt2Rst()

    def get_output_filename(self, path):
        filename, ext = os.path.splitext(path)
        return filename + ".rst"


def main():
    app = Txt2RstConverter()
    app.run()

if __name__ == "__main__":
    main()
