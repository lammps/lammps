/*

  Atz_XML_Parser.cpp
 
  Paul J. Atzberger
  http://atzberger.org/
 
*/

#include "Atz_XML_Parser.h"

const char* Atz_XML_Parser::error_str_code = "Atz_XML_Parser.cpp";

Atz_XML_Parser::Atz_XML_Parser() {
  // TODO Auto-generated constructor stub

}

Atz_XML_Parser::~Atz_XML_Parser() {
  // TODO Auto-generated destructor stub
}


void Atz_XML_Parser::parse(string filename, Atz_XML_SAX_DataHandler *dataHandler) {

  Atz_XML_Parser::parse(filename.c_str(), dataHandler);

}

void Atz_XML_Parser::parse(const char *filename, Atz_XML_SAX_DataHandler *dataHandler) {
	
  const char *error_str_func = "parse()";

  using namespace std;

  ifstream fileStream(filename);
  
  if (fileStream.is_open()) {
    Atz_XML_Parser::parse(&fileStream, dataHandler);
  } else {
  	stringstream message;
    message << "Unable to open the filestream." << endl;
    message << "Filename = " << filename << endl;
    message << "dataHandler = " << dataHandler->DataHandlerName << endl;
    Atz_XML_Package::packageError(error_str_code, error_str_func, message);  
  }

}


void Atz_XML_Parser::parse(ifstream *fileStream, Atz_XML_SAX_DataHandler *dataHandler) {

  using namespace Atz_XML;

  const char *error_str_func = "Atz_XML_Parser()";

  string            line;
  int               count;
  bool              readingTag = false;
  string            myChar;
  string            all;
  string            tag = "";
  stringstream      tagNameStream("");
  stringstream      strStream("");
  string            tagName;
  AttributesType   *tagAttributes = new AttributesType();
  AttributePairType pairStr;
  int               tagType = TAG_TYPE_NULL;
  bool              flagDone = false;
  char              c, c1, c2;

  if (fileStream->is_open()) {

    dataHandler->XML_startDocument(dataHandler);

    count = 0;

    while (fileStream->good()) {

      /* == process any characters before < */
      strStream.str("");
      flagDone = false;
      while (!flagDone) {
        fileStream->get(c);
        if ((!fileStream->good()) || (c == '<')) {
          flagDone = true;
        } else {
          strStream.put(c);
        }
      } /* end while loop */

      /* process data with character callback */
      dataHandler->XML_characters(strStream.str(), dataHandler);

      /* if tag then process using tag parser */
      if (c == '<') {

        /* determine if this is a special tag */
        fileStream->get(c);

        if (c == '!') { /* indicates special tag */
          fileStream->get(c1);
          fileStream->get(c2);
          if ((c1 == '-') && (c2 == '-')) {
            tagType = TAG_TYPE_COMMENT;
          } else {
            fileStream->putback(c1); /* put last char back*/
            fileStream->putback(c2); /* put last char back*/
            tagType = TAG_TYPE_NULL; /* not sure of type (ignore tag) */
          }
          parseUntilEndComment(fileStream); /* parse up to and including symbol '-->'  */

        } else { /* not special tag, so parse as usual */

          /* ==identify the tagName */

          /* put last char back (was used for special checks) */
          fileStream->putback(c);

          /* get the tag name */
          getTagName(fileStream, tagName);

          fileStream->get(c); /* get last char */

          /* if no termination char, then parse attributes */
          tagAttributes->clear();
          if (c == ' ') {
            /* -- collect any tag attributes in hashmap */
            getTagAttributes(fileStream, tagAttributes);
            fileStream->get(c); /* get last char */
          } /* end c == ' ' */

          /* == determine the tag type : start, end, empty */
          tagType = TAG_TYPE_NULL;
          if (c == '>') {
            strStream.str(tagName);
            strStream.get(c);
            if (c == '/') {
              tagType = TAG_TYPE_END;
            } else if (c == '!') {
              strStream.get(c1);
              strStream.get(c2);
              if ((c1 == '-') && (c2 == '-')) {
                tagType = TAG_TYPE_COMMENT;
              }
            } else {
              tagType = TAG_TYPE_START;
            }

          } else if (c == '?') {
            tagType = TAG_TYPE_QUESTMARK;
            fileStream->get(c); /* expect > as next char, otherwise report error */
          } else if (c == '/') {
            tagType = TAG_TYPE_EMPTY;
            fileStream->get(c); /* expect > as next char, otherwise report error */
          } else {
            /* report an error */
          }

          /* == call the tag triggers */
          /* make the tagName passable */
          tagName = getPassableName(tagName); /* clean up modifier characters */

        } /* end of else normal tag */

        /* == trigger action based on the tag type */
        switch (tagType) {

        case TAG_TYPE_START:
          dataHandler->XML_startElement(tagName, tagAttributes, dataHandler);
          break;

        case TAG_TYPE_END:
          dataHandler->XML_endElement(tagName, dataHandler);
          break;

        case TAG_TYPE_EMPTY:
          dataHandler->XML_startElement(tagName, tagAttributes, dataHandler);
          dataHandler->XML_endElement(tagName, dataHandler);
          break;

        case TAG_TYPE_QUESTMARK:
          dataHandler->XML_startElement(tagName, tagAttributes, dataHandler);
          dataHandler->XML_endElement(tagName, dataHandler);
          break;

        case TAG_TYPE_COMMENT:
          /* do nothing */
          break;

        case TAG_TYPE_NULL:
          /* do nothing */
          break;

        default:
          /* report error */
          break;

        } /* end of switch */

      } /* end of tag check */

      count++;

    } /* end of while stream good */


    dataHandler->XML_endDocument(dataHandler);
    fileStream->close();

  } else {
    stringstream message;
    message << "Unable to open the filestream." << endl;
    Atz_XML_Package::packageError(error_str_code, error_str_func, message);
  }

}

void Atz_XML_Parser::removeLeadingWhiteSpace(ifstream*strStream) {
  bool  flagDone = false;
  char  c;

  while (!flagDone) {

    strStream->get(c);

    if ((c == ' ') || (c == '\n')) {
      /* nothing to do, just skip */
    } else {
      flagDone = true;
      strStream->putback(c); /* put non-white back on que */
    }

  } /* end of while loop */

}

void Atz_XML_Parser::getTagName(ifstream *tagStream, string& tagName) {

  bool         flagDone;
  bool         flagFirst = true;
  stringstream tagNameStream("");
  char         c;

  removeLeadingWhiteSpace(tagStream);
  flagDone = false;
  while (!flagDone) {

    tagStream->get(c); /* get next char */

    if ((c == ' ') || (c == '>') || ((c == '/') && (!flagFirst)) ) {
      flagDone = true;
      /* put back on the que */
      tagStream->putback(c);
    } else if (c == '<') {
      /* just skip this char */
    } else {
      tagNameStream.put(c);
    }

    flagFirst = false;

  } /* end while loop */

  tagName = tagNameStream.str(); /* return the tagName */

}


void Atz_XML_Parser::parseUntilEndComment(ifstream *tagStream) {

  bool         flagDone;
  bool         flagFirst = true;
  //stringstream tagNameStream("");
  char         c;
  char         c1;
  char         c2;

  flagDone = false;
  c1 = 0; c2 =0;
  while (!flagDone) {

    tagStream->get(c); /* get next char */

    /* parse until '-->' (assumes no nested comments) */
    if ((c == '>') && (c1 == '-') && (c2 == '-')) {
      flagDone = true;
      /* put back on the que */
      //tagStream->putback(c);
    } else if (c == '<') {
      /* just skip this char */
    } else {
      //tagNameStream.put(c);
    }

    /* record last two characters including current one */
    c2 = c1;
    c1 = c;

    flagFirst = false;

  } /* end while loop */

  //tagName = tagNameStream.str(); /* return the tagName */

}



void   Atz_XML_Parser::getTagAttributes(ifstream *tagStream, Atz_XML::AttributesType *attributes) {

  using namespace Atz_XML;

  bool        flagDone  = false;
  bool        flagValid = true;
  char        c;
  string      attrName;
  string      attrValue;

  /* determine the  name = " " pairs */

  /* read name ignoring white space */
  flagDone = false;
  while (!flagDone) {

    getAttrName(tagStream, attrName); /* handles white space and ";" */

    tagStream->get(c); /* this should be "=", otherwise report error */

    getAttrValue(tagStream, attrValue); /* handles white space and quotes and ";" */

    tagStream->get(c);     /* check next character (should be space or terminal char) */
    tagStream->putback(c); /* put this char back for later routines */

    if ( (!tagStream->good()) || (c == '>') || (c == '/') || (c == '?') ) {
      flagDone = true;
    }

    if (flagValid) {
      /* record name,value pair */
      attributes->insert(AttributePairType(attrName, attrValue));
    }

  } /* end while */

}


void Atz_XML_Parser::getAttrName(ifstream* attrStream, string& attrName) {

  bool         flagDone;
  stringstream attrNameStream("");
  char         c;

  removeLeadingWhiteSpace(attrStream);

  flagDone = false;
  while (!flagDone) {

    attrStream->get(c); /* get next char */

    if (c == '=') {
      flagDone = true;
      /* put back on the que */
      attrStream->putback(c);
    } else if ((c == ' ') || (c == '\n')) {
      /* ignore any spaces */
    } else {
      attrNameStream.put(c);
    }

  } /* end while loop */

  /* return attribute name */
  attrName = attrNameStream.str();

}


void Atz_XML_Parser::getAttrValue(ifstream* attrStream, string& attrValue) {

  bool         flagDone;
  bool         flagRecording = false;
  stringstream attrValueStream("");
  char         c;

  removeLeadingWhiteSpace(attrStream);

  flagDone = false;
  while (!flagDone) {

    attrStream->get(c); /* get next char */

    if (c == '\"') {
      if (!flagRecording) {
        flagRecording = true;
      } else {
        flagRecording = false;
        flagDone = true;
      }
    } else {
      if (flagRecording) {
        attrValueStream.put(c);
      }
    }

  } /* end while loop */

  /* return attribute value */
  attrValue = attrValueStream.str();

}



string Atz_XML_Parser::getPassableName(string tagName) {

  string tagPassedName("");

  /* delete the first character if / */
  if ( (tagName.at(0) == '/') || (tagName.at(0) == '?') ) {
    tagPassedName.assign(tagName, 1, tagName.size() - 1);
  } else {
    tagPassedName.assign(tagName, 0, tagName.size());
  }

  return tagPassedName;
}
