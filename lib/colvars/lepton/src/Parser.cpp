/* -------------------------------------------------------------------------- *
 *                                   Lepton                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the Lepton expression parser originating from              *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2015 Stanford University and the Authors.      *
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

#include "lepton/Parser.h"
#include "lepton/CustomFunction.h"
#include "lepton/Exception.h"
#include "lepton/ExpressionTreeNode.h"
#include "lepton/Operation.h"
#include "lepton/ParsedExpression.h"
#include <cctype>
#include <iostream>

using namespace Lepton;
using namespace std;

static const string Digits = "0123456789";
static const string Operators = "+-*/^";
static const bool LeftAssociative[] = {true, true, true, true, false};
static const int Precedence[] = {0, 0, 1, 1, 3};
static const Operation::Id OperationId[] = {Operation::ADD, Operation::SUBTRACT, Operation::MULTIPLY, Operation::DIVIDE, Operation::POWER};

class Lepton::ParseToken {
public:
    enum Type {Number, Operator, Variable, Function, LeftParen, RightParen, Comma, Whitespace};

    ParseToken(string text, Type type) : text(text), type(type) {
    }
    const string& getText() const {
        return text;
    }
    Type getType() const {
        return type;
    }
private:
    string text;
    Type type;
};

string Parser::trim(const string& expression) {
    // Remove leading and trailing spaces.
    
    int start, end;
    for (start = 0; start < (int) expression.size() && isspace(expression[start]); start++)
        ;
    for (end = (int) expression.size()-1; end > start && isspace(expression[end]); end--)
        ;
    if (start == end && isspace(expression[end]))
        return "";
    return expression.substr(start, end-start+1);
}

ParseToken Parser::getNextToken(const string& expression, int start) {
    char c = expression[start];
    if (c == '(')
        return ParseToken("(", ParseToken::LeftParen);
    if (c == ')')
        return ParseToken(")", ParseToken::RightParen);
    if (c == ',')
        return ParseToken(",", ParseToken::Comma);
    if (Operators.find(c) != string::npos)
        return ParseToken(string(1, c), ParseToken::Operator);
    if (isspace(c)) {
        // White space

        for (int pos = start+1; pos < (int) expression.size(); pos++) {
            if (!isspace(expression[pos]))
                return ParseToken(expression.substr(start, pos-start), ParseToken::Whitespace);
        }
        return ParseToken(expression.substr(start, string::npos), ParseToken::Whitespace);
    }
    if (c == '.' || Digits.find(c) != string::npos) {
        // A number

        bool foundDecimal = (c == '.');
        bool foundExp = false;
        int pos;
        for (pos = start+1; pos < (int) expression.size(); pos++) {
            c = expression[pos];
            if (Digits.find(c) != string::npos)
                continue;
            if (c == '.' && !foundDecimal) {
                foundDecimal = true;
                continue;
            }
            if ((c == 'e' || c == 'E') && !foundExp) {
                foundExp = true;
                if (pos < (int) expression.size()-1 && (expression[pos+1] == '-' || expression[pos+1] == '+'))
                    pos++;
                continue;
            }
            break;
        }
        return ParseToken(expression.substr(start, pos-start), ParseToken::Number);
    }

    // A variable, function, or left parenthesis

    for (int pos = start; pos < (int) expression.size(); pos++) {
        c = expression[pos];
        if (c == '(')
            return ParseToken(expression.substr(start, pos-start+1), ParseToken::Function);
        if (Operators.find(c) != string::npos || c == ',' || c == ')' || isspace(c))
            return ParseToken(expression.substr(start, pos-start), ParseToken::Variable);
    }
    return ParseToken(expression.substr(start, string::npos), ParseToken::Variable);
}

vector<ParseToken> Parser::tokenize(const string& expression) {
    vector<ParseToken> tokens;
    int pos = 0;
    while (pos < (int) expression.size()) {
        ParseToken token = getNextToken(expression, pos);
        if (token.getType() != ParseToken::Whitespace)
            tokens.push_back(token);
        pos += (int) token.getText().size();
    }
    return tokens;
}

ParsedExpression Parser::parse(const string& expression) {
    return parse(expression, map<string, CustomFunction*>());
}

ParsedExpression Parser::parse(const string& expression, const map<string, CustomFunction*>& customFunctions) {
    try {
        // First split the expression into subexpressions.

        string primaryExpression = expression;
        vector<string> subexpressions;
        while (true) {
            string::size_type pos = primaryExpression.find_last_of(';');
            if (pos == string::npos)
                break;
            string sub = trim(primaryExpression.substr(pos+1));
            if (sub.size() > 0)
                subexpressions.push_back(sub);
            primaryExpression = primaryExpression.substr(0, pos);
        }

        // Parse the subexpressions.

        map<string, ExpressionTreeNode> subexpDefs;
        for (int i = 0; i < (int) subexpressions.size(); i++) {
            string::size_type equalsPos = subexpressions[i].find('=');
            if (equalsPos == string::npos)
                throw Exception("subexpression does not specify a name");
            string name = trim(subexpressions[i].substr(0, equalsPos));
            if (name.size() == 0)
                throw Exception("subexpression does not specify a name");
            vector<ParseToken> tokens = tokenize(subexpressions[i].substr(equalsPos+1));
            int pos = 0;
            subexpDefs[name] = parsePrecedence(tokens, pos, customFunctions, subexpDefs, 0);
            if (pos != tokens.size())
                throw Exception("unexpected text at end of subexpression: "+tokens[pos].getText());
        }

        // Now parse the primary expression.

        vector<ParseToken> tokens = tokenize(primaryExpression);
        int pos = 0;
        ExpressionTreeNode result = parsePrecedence(tokens, pos, customFunctions, subexpDefs, 0);
        if (pos != tokens.size())
            throw Exception("unexpected text at end of expression: "+tokens[pos].getText());
        return ParsedExpression(result);
    }
    catch (Exception& ex) {
        throw Exception("Parse error in expression \""+expression+"\": "+ex.what());
    }
}

ExpressionTreeNode Parser::parsePrecedence(const vector<ParseToken>& tokens, int& pos, const map<string, CustomFunction*>& customFunctions,
            const map<string, ExpressionTreeNode>& subexpressionDefs, int precedence) {
    if (pos == tokens.size())
        throw Exception("unexpected end of expression");

    // Parse the next value (number, variable, function, parenthesized expression)

    ParseToken token = tokens[pos];
    ExpressionTreeNode result;
    if (token.getType() == ParseToken::Number) {
        double value;
        stringstream(token.getText()) >> value;
        result = ExpressionTreeNode(new Operation::Constant(value));
        pos++;
    }
    else if (token.getType() == ParseToken::Variable) {
        map<string, ExpressionTreeNode>::const_iterator subexp = subexpressionDefs.find(token.getText());
        if (subexp == subexpressionDefs.end()) {
            Operation* op = new Operation::Variable(token.getText());
            result = ExpressionTreeNode(op);
        }
        else
            result = subexp->second;
        pos++;
    }
    else if (token.getType() == ParseToken::LeftParen) {
        pos++;
        result = parsePrecedence(tokens, pos, customFunctions, subexpressionDefs, 0);
        if (pos == tokens.size() || tokens[pos].getType() != ParseToken::RightParen)
            throw Exception("unbalanced parentheses");
        pos++;
    }
    else if (token.getType() == ParseToken::Function) {
        pos++;
        vector<ExpressionTreeNode> args;
        bool moreArgs;
        do {
            args.push_back(parsePrecedence(tokens, pos, customFunctions, subexpressionDefs, 0));
            moreArgs = (pos < (int) tokens.size() && tokens[pos].getType() == ParseToken::Comma);
            if (moreArgs)
                pos++;
        } while (moreArgs);
        if (pos == tokens.size() || tokens[pos].getType() != ParseToken::RightParen)
            throw Exception("unbalanced parentheses");
        pos++;
        Operation* op = getFunctionOperation(token.getText(), customFunctions);
        try {
            result = ExpressionTreeNode(op, args);
        }
        catch (...) {
            delete op;
            throw;
        }
    }
    else if (token.getType() == ParseToken::Operator && token.getText() == "-") {
        pos++;
        ExpressionTreeNode toNegate = parsePrecedence(tokens, pos, customFunctions, subexpressionDefs, 2);
        result = ExpressionTreeNode(new Operation::Negate(), toNegate);
    }
    else
        throw Exception("unexpected token: "+token.getText());

    // Now deal with the next binary operator.

    while (pos < (int) tokens.size() && tokens[pos].getType() == ParseToken::Operator) {
        token = tokens[pos];
        int opIndex = (int) Operators.find(token.getText());
        int opPrecedence = Precedence[opIndex];
        if (opPrecedence < precedence)
            return result;
        pos++;
        ExpressionTreeNode arg = parsePrecedence(tokens, pos, customFunctions, subexpressionDefs, LeftAssociative[opIndex] ? opPrecedence+1 : opPrecedence);
        Operation* op = getOperatorOperation(token.getText());
        try {
            result = ExpressionTreeNode(op, result, arg);
        }
        catch (...) {
            delete op;
            throw;
        }
    }
    return result;
}

Operation* Parser::getOperatorOperation(const std::string& name) {
    switch (OperationId[Operators.find(name)]) {
        case Operation::ADD:
            return new Operation::Add();
        case Operation::SUBTRACT:
            return new Operation::Subtract();
        case Operation::MULTIPLY:
            return new Operation::Multiply();
        case Operation::DIVIDE:
            return new Operation::Divide();
        case Operation::POWER:
            return new Operation::Power();
        default:
            throw Exception("unknown operator");
    }
}

Operation* Parser::getFunctionOperation(const std::string& name, const map<string, CustomFunction*>& customFunctions) {

    static map<string, Operation::Id> opMap;
    if (opMap.size() == 0) {
        opMap["sqrt"] = Operation::SQRT;
        opMap["exp"] = Operation::EXP;
        opMap["log"] = Operation::LOG;
        opMap["sin"] = Operation::SIN;
        opMap["cos"] = Operation::COS;
        opMap["sec"] = Operation::SEC;
        opMap["csc"] = Operation::CSC;
        opMap["tan"] = Operation::TAN;
        opMap["cot"] = Operation::COT;
        opMap["asin"] = Operation::ASIN;
        opMap["acos"] = Operation::ACOS;
        opMap["atan"] = Operation::ATAN;
        opMap["sinh"] = Operation::SINH;
        opMap["cosh"] = Operation::COSH;
        opMap["tanh"] = Operation::TANH;
        opMap["erf"] = Operation::ERF;
        opMap["erfc"] = Operation::ERFC;
        opMap["step"] = Operation::STEP;
        opMap["delta"] = Operation::DELTA;
        opMap["square"] = Operation::SQUARE;
        opMap["cube"] = Operation::CUBE;
        opMap["recip"] = Operation::RECIPROCAL;
        opMap["min"] = Operation::MIN;
        opMap["max"] = Operation::MAX;
        opMap["abs"] = Operation::ABS;
        opMap["floor"] = Operation::FLOOR;
        opMap["ceil"] = Operation::CEIL;
        opMap["select"] = Operation::SELECT;
    }
    string trimmed = name.substr(0, name.size()-1);

    // First check custom functions.

    map<string, CustomFunction*>::const_iterator custom = customFunctions.find(trimmed);
    if (custom != customFunctions.end())
        return new Operation::Custom(trimmed, custom->second->clone());

    // Now try standard functions.

    map<string, Operation::Id>::const_iterator iter = opMap.find(trimmed);
    if (iter == opMap.end())
        throw Exception("unknown function: "+trimmed);
    switch (iter->second) {
        case Operation::SQRT:
            return new Operation::Sqrt();
        case Operation::EXP:
            return new Operation::Exp();
        case Operation::LOG:
            return new Operation::Log();
        case Operation::SIN:
            return new Operation::Sin();
        case Operation::COS:
            return new Operation::Cos();
        case Operation::SEC:
            return new Operation::Sec();
        case Operation::CSC:
            return new Operation::Csc();
        case Operation::TAN:
            return new Operation::Tan();
        case Operation::COT:
            return new Operation::Cot();
        case Operation::ASIN:
            return new Operation::Asin();
        case Operation::ACOS:
            return new Operation::Acos();
        case Operation::ATAN:
            return new Operation::Atan();
        case Operation::SINH:
            return new Operation::Sinh();
        case Operation::COSH:
            return new Operation::Cosh();
        case Operation::TANH:
            return new Operation::Tanh();
        case Operation::ERF:
            return new Operation::Erf();
        case Operation::ERFC:
            return new Operation::Erfc();
        case Operation::STEP:
            return new Operation::Step();
        case Operation::DELTA:
            return new Operation::Delta();
        case Operation::SQUARE:
            return new Operation::Square();
        case Operation::CUBE:
            return new Operation::Cube();
        case Operation::RECIPROCAL:
            return new Operation::Reciprocal();
        case Operation::MIN:
            return new Operation::Min();
        case Operation::MAX:
            return new Operation::Max();
        case Operation::ABS:
            return new Operation::Abs();
        case Operation::FLOOR:
            return new Operation::Floor();
        case Operation::CEIL:
            return new Operation::Ceil();
        case Operation::SELECT:
            return new Operation::Select();
        default:
            throw Exception("unknown function");
    }
}
