#ifndef LEPTON_PARSER_H_
#define LEPTON_PARSER_H_

/* -------------------------------------------------------------------------- *
 *                                   Lepton                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the Lepton expression parser originating from              *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "windowsIncludes.h"
#include <map>
#include <string>
#include <vector>

namespace Lepton {

class CustomFunction;
class ExpressionTreeNode;
class Operation;
class ParsedExpression;
class ParseToken;

/**
 * This class provides the main interface for parsing expressions.
 */

class LEPTON_EXPORT Parser {
public:
    /**
     * Parse a mathematical expression and return a representation of it as an abstract syntax tree.
     */
    static ParsedExpression parse(const std::string& expression);
    /**
     * Parse a mathematical expression and return a representation of it as an abstract syntax tree.
     *
     * @param customFunctions   a map specifying user defined functions that may appear in the expression.
     *                          The key are function names, and the values are corresponding CustomFunction objects.
     */
    static ParsedExpression parse(const std::string& expression, const std::map<std::string, CustomFunction*>& customFunctions);
private:
    static std::string trim(const std::string& expression);
    static std::vector<ParseToken> tokenize(const std::string& expression);
    static ParseToken getNextToken(const std::string& expression, int start);
    static ExpressionTreeNode parsePrecedence(const std::vector<ParseToken>& tokens, int& pos, const std::map<std::string, CustomFunction*>& customFunctions,
            const std::map<std::string, ExpressionTreeNode>& subexpressionDefs, int precedence);
    static Operation* getOperatorOperation(const std::string& name);
    static Operation* getFunctionOperation(const std::string& name, const std::map<std::string, CustomFunction*>& customFunctions);
};

} // namespace Lepton

#endif /*LEPTON_PARSER_H_*/
