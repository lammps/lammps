txt2html - a text to HTML conversion tool :h3

[txt2html] is a simple tool for converting text files into HTML files.
Text files can contain simple formatting and mark-up commands that
[txt2html] converts into HTML.

[txt2html] was written by "Steve Plimpton"_sjp.  I use it for
"documentation"_doc and "WWW pages"_www.  Anna Reese added the table
formatting options.

See the "example.txt"_example.txt and "example.html"_example.html
files in the [txt2html] directory for examples of what all the
formatting commands and mark-up syntax end up looking like in HTML.

:link(sjp,http://www.cs.sandia.gov/~sjplimp)
:link(www,http://www.cs.sandia.gov/~sjplimp)
:link(doc,http://www.cs.sandia.gov/~sjplimp/lammps.html)

:line

[Syntax:]

txt2html file
  read from text file, write HTML to standard output
txt2html file1 file2 file3 ...
  read each argument as text file, write one HTML file per argument :dl

Input files are first opened with the specified name.  If that fails,
a ".txt" suffix is added.  Output files are created with an ".html"
suffix, which is either added or replaces the ".txt" suffix.

:line

[Compiling:]

The source for [txt2html] is a single C++ file.  Compile it by typing:

g++ -o txt2html txt2html.cpp :pre

:line

[How the tool works:]

[txt2html] reads a text file, one {paragraph} at a time.  A paragraph
ends with:

  a blank line
  a line whose final word starts with ":" (a format string)
  the end of the file :ul

Any line in the paragraph which ends with "\" is concatenated to the
following line by removing the "\" character and following newline.
This can be useful for some of the formatting commands described below
that operate on individual lines in the paragraph.

If a paragraph starts with a "&lt;" character and ends with a "&gt;"
character, it is treated as raw HTML and is written directly into the
output file.

If a paragraph does not end with a format string, then it is
surrounded with HTML paragraph markers (&lt;P&gt; and &lt;/P&gt;),
"mark-up"_#markup is performed, and the paragraph is written to the
output file.

If the paragraph ends with a format string, then "formatting"_#format
is performed, "mark-up"_#markup is performed, and the paragraph is
written to the output file.

:line

[Formatting:] :link(format)

A format string is the last word of a paragraph if it starts with a
":" character.  A format string contains one or more comma-separated
commands, like ":ulb,l" or ":c,h3".  Note that a format string cannot
contain spaces, else it would not be the last word.  An individual
command can have 0 or more arguments:

  {b} or {line()} = 0 arguments
  {image(file)} = 1 argument
  {link(alias,value)} = 2 or more comma-separated arguments :ul

Format commands add HTML markers at the beginning or end of the
paragraph and individual lines.  Commands are processed in the order
they appear in the format string.  Thus if two commands add HTML
markers to the beginning of the paragraph, the 2nd command's marker
will appear 2nd.  The reverse is true at the end of the paragraph; the
2nd command's marker will appear 1st.  Some comands, like {line} or
{image} make most sense if used as stand-alone commands without an
accompanying paragraph.

Commands that format the entire paragraph:

  p --&gt; surround the paragraph with &lt;P&gt; &lt;/P&gt;
  b --&gt; put &lt;BR&gt; at the end of the paragraph
  pre --&gt; surround the paragraph with &lt;PRE&gt; &lt;/PRE&gt;
  c --&gt; surround the paragraph with &lt;CENTER&gt; &lt;/CENTER&gt;
  h1,h2,h3,h4,h5,h6 --&gt; surround the paragraph with \
                           &lt;H1&gt; &lt;/H1&gt;, etc :ul

Commands that format the lines of the paragraph as a list:

  ul --&gt; surround the paragraph with &lt;UL&gt; &lt;/UL&gt;, \
    put &lt;LI&gt; at start of every line
  ol --&gt; surround the paragraph with &lt;OL&gt; &lt;/OL&gt;, \
    put &lt;LI&gt; at start of every line
  dl --&gt; surround the paragraph with &lt;DL&gt; &lt;/DL&gt;, \
    alternate &lt;DT&gt; and &lt;DD&gt; at start of every line :ul

Commands that treat the paragraph as one entry in a list:

  l --&gt; put &lt;LI&gt; at the beginning of the paragraph
  dt --&gt; put &lt;DT&gt; at the beginning of the paragraph
  dd --&gt; put &lt;DD&gt; at the beginning of the paragraph
  ulb --&gt; put &lt;UL&gt; at the beginning of the paragraph
  ule --&gt; put &lt;/UL&gt; at the end of the paragraph
  olb --&gt; put &lt;OL&gt; at the beginning of the paragraph
  ole --&gt; put &lt;/OL&gt; at the end of the paragraph
  dlb --&gt; put &lt;DL&gt; at the beginning of the paragraph
  dle --&gt; put &lt;/DL&gt; at the end of the paragraph :ul

Commands applied to each line of the paragraph:

  all(p) --&gt; surround each line with &lt;P&gt; &lt;/P&gt;
  all(c) --&gt; surround each line with &lt;CENTER&gt; &lt;/CENTER&gt;
  all(b) --&gt; append a &lt;BR&gt; to each line
  all(l) --&gt; prepend a &lt;LI&gt; to each line :ul

Special commands (all HTML is inserted at beginning of paragraph):

  line --&gt; insert a horizontal line = &lt;HR&gt;
  image(file) --&gt; insert an image = &lt;IMG SRC = "file"&gt;
  image(file,link) --&gt; insert an image that when clicked on goes to link
  link(name) --&gt; insert a named link that can be referred to \
    elsewhere (see "mark-up"_#markup) = &lt;A NAME = "name"&gt;&lt;/A&gt;
  link(alias,value) --&gt; define a link alias that can be used \
    elsewhere in this file (see "mark-up"_#markup) :ul


Table command:

  tb(c=3,b=5,w=100%,a=c) --&gt; format the paragraph as a table :ul

Arguments within tb() can appear in any order and are all optional,
since they each have default values.

  c=N --&gt; Make an N-column table.  Treat the paragraph as one
  long list of entries (separated by the separator character) and put
  them into N columns one after the other.  If N = 0, treat each line
  of the paragraph as one row of the table with as many columns as
  there are maximum entries in any line.  Default is c=0. :ulb,l

  s=: --&gt; Use the character string following the equal sign as
  the separator between entries.  Default separator is a comma "," which
  you cannot specify directly since the comma delimits the tb() arguments :l

  b=N --&gt; Create a border N pixels wide.  If N is 0, there is no
  border between or outside the cells.  If N is 1, there is a minimal
  border between and outside all cells.  For N > 1, the border between
  cells does not change but the outside border gets wider.  Default is
  b=1. :l

  w=N or w=N% --&gt The first form makes each cell of the table at
  least N pixels wide.  The second form makes the entire table take up
  N% of the width of the browser window.  Default is w=0 which means
  each cell will be just as wide as the text it contains.  :l

  a=X --&gt Align the entire table at the left, center, or right of the
  browser window, for X = "l", "c", or "r".  Default is a=c. :l

  ea=X --&gt Align the text in each entry at the left, center, or
  right of its cell, for X = "l", "c", or "r".  Default is browser's 
  default (typically left). :l

  eva=X --&gt Vertically align the text in each entry at the
  top, middle, baseline, or bottom of its cell, for X = "t", "m", "ba", 
  or "bo".  Default is browser's default (typically middle). :l

  cwM=N or cwM=N% --&gt The first form makes column M be at least
  N pixels wide.  The second form makes column M take up N% of the
  width of the browser window.  This setting overrides the "w"
  argument for column M.  Only one column per table can be tweaked
  with this argument.  Default is no settings for any column. :l

  caM=X --&gt Align the text in each entry of column M at the left,
  center, or right of its cell, for X = "l", "c", or "r".  This
  setting overrides the "ea" argument for column M.  Only one column
  per table can be tweaked with this argument.  Default is no settings
  for any column. :l

  cvaM=X --&gt Vertically align the text in each entry of column m
  at the top, middle, baseline, or bottom of its cell, for X = "t",
  "m", "ba", or "bo".  This setting overrides the "eva" argument for
  column M.  Only one column per table can be tweaked with this
  argument.  Default is no settings for any column. :l,ule

:line

[Mark-up:] :link(markup)

The text of the paragraph is scanned for special mark-up characters
which are converted into HTML.

Bold and italic characters:

  <UL> <LI> "[" (left brace) --&gt; turn-on bold by inserting a &lt;B&gt;
  <LI> "]" (right brace) --&gt; turn-off bold by inserting a &lt;/B&gt;
  <LI> "{" (left bracket) --&gt; turn-on italics by inserting a &lt;I&gt;
  <LI> "}" (right bracket) --&gt; turn-off italics by 
    inserting a &lt;/I&gt; </UL>

If a backslash ESCAPE character '\' preceeds any of the bold/italic
mark-up characters, then mark-up is not performed; the mark-up
character is simply left in the text.

Links are inserted by enclosing a section of text in double quotes,
and appending an underscore to the ending quote, followed by the link.
The link ends when whitespace is found, except that trailing
punctuation characters (comma, period, semi-colon, colon, question
mark, exclamation point, parenthesis) are not considered part of the
link.

<P> A link of the form "text"_link becomes &lt;A HREF =
"link"&gt;text&lt;/A&gt; in the HTML output.  The only exception is if
"link" is defined elsewhere in the file as an alias (see the link
command above).  In that case, the value is used instead of the alias
name. </P>

With these rules, links can take several forms.

<UL> <LI> "This links"_#abc to another part of this file which is
labeled with a :link(abc) command. <BR>
<LI> "This links"_other.html to another file named other.html. <BR>
<LI> "This links"_other.html#abc to another file which has an "abc"
location defined internally. <BR>
<LI> "This links"_http://www.google.com to a WWW site. <BR>
<LI> "This"_M12 could be used in place of any of the above forms.  It
requires an alias like :link(M12,http://www.google.com) to be defined
elsewhere in the file. </UL>
