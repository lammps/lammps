"""Contains the functions used to read the input file and print the checkpoint
files with xml formatting.

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


Functions:
   xml_node: Class to handle a particular xml tag.
   xml_handler: Class giving general xml data reading methods. 
   xml_parse_string: Parses a string made from a section of a xml input file.
   xml_parse_file: Parses an entire xml input file.
   read_type: Reads a string and outputs data of a specified type.
   read_float: Reads a string and outputs a float.
   read_int: Reads a string and outputs an integer.
   read_bool: Reads a string and outputs a boolean.
   read_list: Reads a string and outputs a list.
   read_array: Reads a string and outputs an array.
   read_tuple: Reads a string and outputs a tuple.
   read_dict: Reads a string and outputs a dictionary.
   write_type: Writes a string from data of a specified type.
   write_list: Writes a string from a list.
   write_tuple: Writes a string from a tuple.
   write_float: Writes a string from a float.
   write_bool: Writes a string from a boolean.
   write_dict: Writes a string from a dictionary.
"""

__all__ = ['xml_node', 'xml_handler', 'xml_parse_string', 'xml_parse_file',
           'read_type', 'read_float', 'read_int', 'read_bool', 'read_list',
           'read_array', 'read_tuple', 'read_dict', 'write_type', 'write_list',
           'write_tuple', 'write_float', 'write_bool', 'write_dict']

from xml.sax import parseString, parse
from xml.sax.handler import ContentHandler 
import numpy as np
import string

class xml_node(object):
   """Class to handle a particular xml tag.

   Tags are generally written in the form 
   <tag_name attribs="attrib_data"> main_data </tag_name>. This class holds
   tag_name, attrib_data and main_data separately so they can be used to
   create the objects with the appropriate names and data.

   Attributes:
      attribs: The attribute data for the tag.
      fields: The rest of the data.
      name: The tag name. 
   """

   def __init__(self, attribs=None, name="", fields=None):
      """Initialises xml_node.

      Args:
         attribs: An optional dictionary giving attribute data. Defaults to {}.
         fields: An optional dictionary holding all the data between the start
            and end tags, including information about other nodes. 
            Defaults to {}.
         name: An optional string giving the tag name. Defaults to ''.
      """

      if attribs is None:
         attribs = {}
      if fields is None:
         fields = []

      self.attribs = attribs
      self.name = name         
      self.fields = fields


class xml_handler(ContentHandler):
   """Class giving general xml_reading methods.

   Uses the standard python xml_reader to read the different kinds of data.
   Keeps track of the heirarchial nature of an xml file by recording the level
   of nesting, so that the correct data and attributes can be associated with 
   the correct tag name. 

   Attributes:
      root: An xml_node object for the root node.
      open: The list of the tags that the parser is currently between the start
         and end tags of.
      level: The level of nesting that the parser is currently at.
      buffer: A list of the data found between the tags at the different levels
         of nesting.
   """

   def __init__(self):
      """Initialises xml_handler."""

      #root xml node with all the data
      self.root = xml_node(name="root", fields=[])
      self.open = [self.root]
      #current level of the hierarchy
      self.level = 0
      #Holds all the data between each of the tags.
      #If level = 1, then buffer[0] holds all the data collected between the 
      #root tags, and buffer[1] holds all the data collected between the 
      #first child tag.
      self.buffer = [[""]]
      
   def startElement(self, name, attrs): 
      """Reads an opening tag.

      Adds the opening tag to the list of open tags, adds a new space in the
      buffer, reads the appropriate attributes and adds a new level to the 
      heirarchy.

      Args:
         name: The tag_name.
         attrs: The attribute data.
      """

      #creates a new node
      newnode = xml_node(attribs=dict((k,attrs[k]) for k in attrs.keys()), name=name, fields=[])
      #adds it to the list of open nodes
      self.open.append(newnode)
      #adds it to the list of fields of the parent tag
      self.open[self.level].fields.append((name,newnode))
      #gets ready to read new data
      self.buffer.append([""])
      self.level += 1      

   def characters(self, data):
      """Reads data.

      Adds the data to the buffer of the current level of the heirarchy.
      Data is read as a string, and needs to be converted to the required
      type later.

      Args:
         data: The data to be read.
      """

      self.buffer[self.level].append(data)

   def endElement(self, name):
      """Reads a closing tag.

      Once all the data has been read, and the closing tag found, the buffer
      is read into the appropriate field.

      Args:
         name: The tag_name.
      """

      #all the text found between the tags stored in the appropriate xml_node
      #object
      self.buffer[self.level] = ''.join(self.buffer[self.level])
      self.open[self.level].fields.append(("_text" , self.buffer[self.level]))
      #'closes' the xml_node object, as we are no longer within its tags, so
      #there is no more data to be added to it.
      #Note that the xml_node is still held within the parent tag, so we
      #no longer require this xml node object.
      self.buffer.pop(self.level)
      self.open.pop(self.level)
      self.level -= 1

def xml_parse_string(buf):
   """Parses a string made from a section of a xml input file.

   Args:
      buf: A string in correct xml format.

   Returns:
      A xml_node for the root node of the file.
   """

   myhandle = xml_handler()
   parseString(buf, myhandle)
   return myhandle.root

def xml_parse_file(stream):
   """Parses an entire xml input file.

   Args:
      stream: A string describing a xml formatted file.

   Returns:
      A xml_node for the root node of the file.
   """

   myhandle = xml_handler()
   parse(stream, myhandle)
   return myhandle.root

def read_type(type, data):
   """Reads a string and outputs data of a specified type.

   Args:
      type: The data type of the target container.
      data: The string to be read in.

   Raises:
      TypeError: Raised if it tries to read into a data type that has not been
         implemented.

   Returns:
      An object of type type.
   """

   if not type in readtype_funcs:
      raise TypeError("Conversion not available for given type")
   return type(readtype_funcs[type](data))

def read_float(data):
   """Reads a string and outputs a float.

   Args:
      data: The string to be read in.

   Raises: 
      ValueError: Raised if the input data is not of the correct format.

   Returns:
      A float.
   """

   return float(data)

def read_int(data):
   """Reads a string and outputs a integer.

   Args:
      data: The string to be read in.

   Raises: 
      ValueError: Raised if the input data is not of the correct format.

   Returns:
      An integer.
   """

   return int(data)

def read_bool(data):
   """Reads a string and outputs a boolean.

   Takes a string of the form 'true' or 'false', and returns the appropriate
   boolean.

   Args:
      data: The string to be read in.

   Raises:
      ValueError: Raised if the string is not 'true' or 'false'.

   Returns:
      A boolean.
   """


   if data.strip().upper() == "TRUE":
      return True
   elif data.strip().upper() == "FALSE":
      return False
   else:
      raise ValueError(data + " does not represent a bool value")

def read_list(data, delims="[]", split=",", strip=" \n\t'"):
   """Reads a formatted string and outputs a list.

   The string must be formatted in the correct way.
   The start character must be delimiters[0], the end character
   must be delimiters[1] and each element must be split along 
   the character split. Characters at the beginning or
   end of each element in strip are ignored. The standard list format is of the
   form '[array[0], array[1],..., array[n]]', which is used for actual lists.
   Other formats are used for tuples and dictionaries.

   Args:
      data: The string to be read in. '[]' by default.
      delims: A string of two characters giving the first and last character of
         the list format. ',' by default.
      split: The character between different elements of the list format.
      strip: Characters to be removed from the beginning and end of each 
         element. ' \n\t' by default.

   Raises: 
      ValueError: Raised if the input data is not of the correct format.

   Returns:
      A list of strings.
   """

   try:
      begin = data.index(delims[0])
      end = data.index(delims[1])
   except ValueError:
      raise ValueError("Error in list syntax: could not locate delimiters")
   
   rlist = data[begin+1:end].split(split)
   for i in range(len(rlist)):
      rlist[i] = rlist[i].strip(strip)

   # handles empty lists correctly
   if len(rlist) == 1 and rlist[0] == "":
      rlist = [] 
    
   return rlist

def read_array(dtype, data):
   """Reads a formatted string and outputs an array.

   The format is as for standard python arrays, which is
   [array[0], array[1], ... , array[n]]. Note the use of comma separators, and 
   the use of square brackets.

   Args:
      data: The string to be read in.
      dtype: The data type of the elements of the target array.

   Raises: 
      ValueError: Raised if the input data is not of the correct format.

   Returns:
      An array of data type dtype.
   """

   rlist = read_list(data)
   for i in range(len(rlist)):
      rlist[i] = read_type(dtype,rlist[i])
   
   return np.array(rlist, dtype)

def read_tuple(data, delims="()", split=",", strip=" \n\t'", arg_type=int):
   """Reads a formatted string and outputs a tuple.

   The format is as for standard python tuples, which is
   (tuple[0], tuple[1], ... , tuple[n]). Note the comma 
   separators, and the use of brackets.

   Args:
      data: The string to be read in.
      delims: A string of two characters giving the first and last character of
         the list format. ',' by default.
      split: The character between different elements of the list format.
      strip: Characters to be removed from the beginning and end of each 
         element. ' \n\t' by default.
      arg_type: The strings in the input will be converted, and a tuple 
         of ar_type will be returned. 

   Raises: 
      ValueError: Raised if the input data is not of the correct format.

   Returns:
      A tuple of elements of the specified data type.
   """

   rlist = read_list(data, delims=delims, split=split, strip=strip)
   return tuple([arg_type(i) for i in rlist])

def read_dict(data, delims="{}", split=",", key_split=":", strip=" \n\t"):
   """Reads a formatted string and outputs a dictionary.

   The format is as for standard python dictionaries, which is
   {keyword[0]: arg[0], keyword[1]: arg[1], ... , keyword[n]: arg[n]}. Note the
   comma separators, and the use of curly brackets.

   Args:
      data: The string to be read in.
      delims: A string of two characters giving the first and last character of
         the list format. ',' by default.
      split: The character between different elements of the list format.
      key_split: The character between the key word and the value.
      strip: Characters to be removed from the beginning and end of each 
         element. ' \n\t' by default.

   Raises: 
      ValueError: Raised if the input data is not of the correct format.

   Returns:
      A dictionary of strings.
   """

   rlist = read_list(data, delims=delims, split=split, strip=strip)
   def mystrip(data):
      return data.strip(strip)
   rdict = {}
   for s in rlist:
      rtuple = map(mystrip,s.split(key_split))      
      if not len(rtuple) == 2:
         raise ValueError("Format for a key:value format is wrong for item " + s)
      rdict[rtuple[0]] = rtuple[1]
      
   return rdict   
      
readtype_funcs = {np.ndarray: read_array, dict: read_dict, float: read_float, int: read_int, bool: read_bool, str: string.strip, tuple: read_tuple, np.uint : read_int}

def write_type(type, data):
   """Writes a formatted string from a value of a specified type.

   Args:
      type: The data type of the value.
      data: The value to be read in.

   Raises:
      TypeError: Raised if it tries to write from a data type that has not been
         implemented.

   Returns:
      A formatted string.
   """

   if not type in writetype_funcs:
      raise TypeError("Conversion not available for given type")
   return writetype_funcs[type](data)

def write_list(data, delims="[]"):
   """Writes a formatted string from a list.

   The format of the output is as for a standard python list, 
   [list[0], list[1],..., list[n]]. Note the space after the commas, and the
   use of square brackets.

   Args:
      data: The value to be read in.
      delims: An optional string of two characters giving the first and last
         character to be printed. Defaults to "[]".

   Returns:
      A formatted string.
   """

   rstr = delims[0]
   
   for v in data:
      rstr += str(v) + ", "
   
   rstr = rstr.rstrip(", ")
   rstr += delims[1]
   return rstr

def write_tuple(data):
   """Writes a formatted string from a tuple.

   The format of the output is as for a standard python tuple, 
   (tuple[0], tuple[1],..., tuple[n]). Note the space after the commas, and the
   use of brackets.

   Args:
      data: The value to be read in.

   Returns:
      A formatted string.
   """

   return write_list(data, delims="()")

def write_float(data):
   """Writes a formatted string from a float.

   Floats are printed out in exponential format, to 8 decimal places and 
   filling up any spaces under 16 not used with spaces.

   For example 1.0 --> '  1.00000000e+00'

   Args:
      data: The value to be read in.

   Returns:
      A formatted string.
   """

   return "%16.8e" % (data)

def write_bool(data):
   """Writes a formatted string from a float.

   Booleans are printed as a string of either ' true' or 'false'. Note that
   both are printed out as exactly 5 characters.

   Args:
      data: The value to be read in.

   Returns:
      A formatted string.
   """

   return "%5.5s" % (str(data))

def write_dict(data, delims="{}"):
   """Writes a formatted string from a dictionary.

   The format of the output is as for a standard python dictionary, 
   {keyword[0]: arg[0], keyword[1]: arg[1],..., keyword[n]: arg[n]}. Note the 
   space after the commas, and the use of curly brackets.

   Args:
      data: The value to be read in.
      delims: An optional string of two characters giving the first and last
         character to be printed. Defaults to "{}".

   Returns:
      A formatted string.
   """

   rstr = delims[0]
   for v in data:
      rstr += str(v) + ": " + str(data[v]) + ", "
   rstr = rstr.strip(", ")
   rstr += delims[1]
   return rstr

writetype_funcs = {float: write_float, dict: write_dict, int: str, bool: write_bool, str: string.strip, tuple: write_tuple, np.uint : str}
