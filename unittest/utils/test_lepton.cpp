
// Adapted for GoogleTest from TestParser.cpp from OpenMM

#include "lammps.h"

#include "info.h"
#include "input.h"
#include "update.h"
#include "variable.h"

#include "../../src/LEPTON/lepton_utils.h"
#include "Lepton.h"
#include "lepton/CompiledVectorExpression.h"
#include "utils.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "../testing/core.h"

#include <exception>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>

using LAMMPS_NS::utils::split_words;
using ::testing::StrEq;

bool verbose = false;

class LeptonUtilsTest : public LAMMPSTest {
protected:
    LAMMPS_NS::Variable *variable;

    void SetUp() override
    {
        testbinary = "LeptonUtilsTest";
        args       = {"-log", "none", "-echo", "screen", "-nocite", "-v", "num", "1"};
        LAMMPSTest::SetUp();
        BEGIN_HIDE_OUTPUT();
        command("region box block 0 1 0 1 0 1");
        command("create_box 1 box");
        END_HIDE_OUTPUT();
        variable = lmp->input->variable;
    }
};

// remove quotes and spaces from expression

TEST(LeptonUtils, condense)
{
    ASSERT_THAT(LeptonUtils::condense("\"k*r^2; k=250.0\""), StrEq("k*r^2;k=250.0"));
    ASSERT_THAT(LeptonUtils::condense("'k2*r^2 + k3*r^3 + k4*r^4; k2=300.0; k3=-100.0; k4=50.0'"),
                StrEq("k2*r^2+k3*r^3+k4*r^4;k2=300.0;k3=-100.0;k4=50.0"));
    ASSERT_THAT(LeptonUtils::condense("k*(r-0.2)^2;k=500.0"), StrEq("k*(r-0.2)^2;k=500.0"));
    ASSERT_THAT(LeptonUtils::condense("\"xx' \"'xx"), StrEq("xxxx"));
    ASSERT_THAT(LeptonUtils::condense("\t \"x\n\r"), StrEq("x"));
}

// substitute variable references (v_<name>) with values

TEST_F(LeptonUtilsTest, substitute)
{
    BEGIN_HIDE_OUTPUT();
    command("variable val1 index 100.0");
    command("variable pre equal 0.001*step");
    END_HIDE_OUTPUT();
    ASSERT_THAT(LeptonUtils::substitute("v_num", lmp), StrEq("1"));
    ASSERT_THAT(LeptonUtils::substitute("eps*v_val1*k", lmp), StrEq("eps*100.0*k"));
    ASSERT_THAT(LeptonUtils::substitute("(2.5/v_pre)", lmp), StrEq("(2.5/0)"));
    lmp->update->reset_timestep(100LL, false);
    ASSERT_THAT(LeptonUtils::substitute("(2.5/v_pre)", lmp), StrEq("(2.5/0.1)"));

    bool caught = false;
    try {
        LeptonUtils::substitute("v_none", lmp);
    } catch (std::exception &e) {
        ASSERT_THAT(e.what(), StrEq("Variable none in expression v_none does not exist"));
        caught = true;
    }
    ASSERT_TRUE(caught);
}

// zbl() custom function

TEST(LeptonCustomFunction, zbl)
{
    Lepton::ZBLFunction zbl(1.0, 1.0, 1.0);
    std::map<std::string, Lepton::CustomFunction *> functions = {std::make_pair("zbl", &zbl)};
    std::map<std::string, double> variables = {std::make_pair("zi", 6), std::make_pair("zj", 6),
                                               std::make_pair("r", 2.0)};

    auto parsed = Lepton::Parser::parse("zbl(zi, zj, r)", functions);
    auto zbldzi = parsed.differentiate("zi");
    auto zbldzj = parsed.differentiate("zj");
    auto zbldr  = parsed.differentiate("r");
    auto zbld2r = zbldr.differentiate("r");

    double value = parsed.evaluate(variables);
    ASSERT_DOUBLE_EQ(value, 0.065721538245489763);
    value = zbldr.evaluate(variables);
    ASSERT_DOUBLE_EQ(value, -0.15481915325334394);
    variables["r"] = 1.0;
    value          = parsed.evaluate(variables);
    ASSERT_DOUBLE_EQ(value, 1.0701488641432269);
    value = zbldr.evaluate(variables);
    ASSERT_DOUBLE_EQ(value, -3.6376386525054412);
    variables["zi"] = 13.0;
    value           = parsed.evaluate(variables);
    ASSERT_DOUBLE_EQ(value, 1.8430432789454971);
    value = zbldr.evaluate(variables);
    ASSERT_DOUBLE_EQ(value, -6.5373118484557642);
    variables["zj"] = 13.0;
    value           = parsed.evaluate(variables);
    ASSERT_DOUBLE_EQ(value, 3.1965196467438446);
    value = zbldr.evaluate(variables);
    ASSERT_DOUBLE_EQ(value, -11.804490148948526);

    // check for unsupported derivatives
    ASSERT_ANY_THROW(value = zbldzi.evaluate(variables));
    ASSERT_ANY_THROW(value = zbldzj.evaluate(variables));
    ASSERT_ANY_THROW(value = zbld2r.evaluate(variables));
}

/**
 * This is a custom function equal to f(x,y) = 2*x*y.
 */

class ExampleFunction : public Lepton::CustomFunction {
    int getNumArguments() const { return 2; }
    double evaluate(const double *arguments) const { return 2.0 * arguments[0] * arguments[1]; }
    double evaluateDerivative(const double *arguments, const int *derivOrder) const
    {
        if (derivOrder[0] == 1) {
            if (derivOrder[1] == 0)
                return 2.0 * arguments[1];
            else if (derivOrder[1] == 1)
                return 2.0;
        }
        if (derivOrder[1] == 1 && derivOrder[0] == 0) return 2.0 * arguments[0];
        return 0.0;
    }
    Lepton::CustomFunction *clone() const { return new ExampleFunction(); }
};

/**
 * Verify that an expression gives the correct value.
 */

void verifyEvaluation(const std::string &expression, double expectedValue)
{
    std::map<std::string, Lepton::CustomFunction *> customFunctions;
    Lepton::ParsedExpression parsed = Lepton::Parser::parse(expression, customFunctions);
    double value                    = parsed.evaluate();
    ASSERT_NEAR(expectedValue, value, 1e-10);

    // Try optimizing it and make sure the result is still correct.

    value = parsed.optimize().evaluate();
    ASSERT_NEAR(expectedValue, value, 1e-10);

    // Create an ExpressionProgram and see if that also gives the same result.

    Lepton::ExpressionProgram program = parsed.createProgram();
    value                             = program.evaluate();
    ASSERT_NEAR(expectedValue, value, 1e-10);

    // Create a CompiledExpression and see if that also gives the same result.

    Lepton::CompiledExpression compiled = parsed.createCompiledExpression();
    value                               = compiled.evaluate();
    ASSERT_NEAR(expectedValue, value, 1e-10);
}

/**
 * Verify that an expression with variables gives the correct value.
 */

void verifyEvaluation(const std::string &expression, double x, double y, double expectedValue)
{
    if (verbose) std::cout << "Checking expression: " << expression << "\n";
    std::map<std::string, double> variables;
    variables["x"]                  = x;
    variables["y"]                  = y;
    Lepton::ParsedExpression parsed = Lepton::Parser::parse(expression);
    double value                    = parsed.evaluate(variables);
    ASSERT_NEAR(expectedValue, value, 1e-10);

    // Try optimizing it and make sure the result is still correct.

    value = parsed.optimize().evaluate(variables);
    ASSERT_NEAR(expectedValue, value, 1e-10);

    // Try optimizing with predefined values for the variables.

    value = parsed.optimize(variables).evaluate();
    ASSERT_NEAR(expectedValue, value, 1e-10);

    // Create an ExpressionProgram and see if that also gives the same result.

    Lepton::ExpressionProgram program = parsed.createProgram();
    value                             = program.evaluate(variables);
    ASSERT_NEAR(expectedValue, value, 1e-10);

    // Create a CompiledExpression and see if that also gives the same result.

    Lepton::CompiledExpression compiled = parsed.createCompiledExpression();
    if (compiled.getVariables().find("x") != compiled.getVariables().end())
        compiled.getVariableReference("x") = x;
    if (compiled.getVariables().find("y") != compiled.getVariables().end())
        compiled.getVariableReference("y") = y;
    value = compiled.evaluate();
    ASSERT_NEAR(expectedValue, value, 1e-10);

    // Try specifying memory locations for the compiled expression.

    std::map<std::string, double *> variablePointers;
    variablePointers["x"]                = &x;
    variablePointers["y"]                = &y;
    Lepton::CompiledExpression compiled2 = parsed.createCompiledExpression();
    compiled2.setVariableLocations(variablePointers);
    value = compiled2.evaluate();
    ASSERT_NEAR(expectedValue, value, 1e-10);
    ASSERT_EQ(&x, &compiled2.getVariableReference("x"));
    ASSERT_EQ(&y, &compiled2.getVariableReference("y"));

    // Try evaluating it as a vector.

    for (int width : Lepton::CompiledVectorExpression::getAllowedWidths()) {
        Lepton::CompiledVectorExpression vector = parsed.createCompiledVectorExpression(width);
        for (int i = 0; i < width; i++) {
            if (vector.getVariables().find("x") != vector.getVariables().end())
                for (int j = 0; j < width; j++)
                    vector.getVariablePointer("x")[j] = (i == j ? x : -100.0);
            if (vector.getVariables().find("y") != vector.getVariables().end())
                for (int j = 0; j < width; j++)
                    vector.getVariablePointer("y")[j] = (i == j ? y : -100.0);
            const float *result = vector.evaluate();
            ASSERT_NEAR(expectedValue, result[i], 1e-6);
        }
    }

    // Specify memory locations for the vector expression.

    float xvec[8], yvec[8];
    std::map<std::string, float *> vecVariablePointers;
    vecVariablePointers["x"] = xvec;
    vecVariablePointers["y"] = yvec;
    for (int width : Lepton::CompiledVectorExpression::getAllowedWidths()) {
        Lepton::CompiledVectorExpression vector2 = parsed.createCompiledVectorExpression(width);
        vector2.setVariableLocations(vecVariablePointers);
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < width; j++) {
                xvec[j] = (i == j ? x : -100.0);
                yvec[j] = (i == j ? y : -100.0);
            }
            const float *result = vector2.evaluate();
            ASSERT_NEAR(expectedValue, result[i], 1e-6);
        }
    }

    // Make sure that variable renaming works.

    variables.clear();
    variables["w"] = x;
    variables["y"] = y;
    std::map<std::string, std::string> replacements;
    replacements["x"] = "w";
    value             = parsed.renameVariables(replacements).evaluate(variables);
    ASSERT_NEAR(expectedValue, value, 1e-10);
}

/**
 * Confirm that a parse error gets thrown.
 */

void verifyInvalidExpression(const std::string &expression)
{
    if (verbose) std::cout << "Checking invalid expression: " << expression << "\n";
    try {
        Lepton::Parser::parse(expression);
    } catch (const std::exception &) {
        return;
    }
    throw std::exception();
}

/**
 * Verify that two numbers have the same value.
 */

void assertNumbersEqual(double val1, double val2, double tol = 1e-10)
{
    const double inf = std::numeric_limits<double>::infinity();
    if (val1 == val1 || val2 == val2)           // If both are NaN, that's fine.
        if (val1 != inf || val2 != inf)         // Both infinity is also fine.
            if (val1 != -inf || val2 != -inf) { // Same for -infinity.
                ASSERT_NEAR(val1, val2, tol);
            }
}

/**
 * Verify that two expressions give the same value.
 */

void verifySameValue(const Lepton::ParsedExpression &exp1, const Lepton::ParsedExpression &exp2,
                     double x, double y)
{
    std::map<std::string, double> variables;
    variables["x"] = x;
    variables["y"] = y;
    double val1    = exp1.evaluate(variables);
    double val2    = exp2.evaluate(variables);
    assertNumbersEqual(val1, val2);

    // Now create CompiledExpressions from them and see if those also match.

    Lepton::CompiledExpression compiled1 = exp1.createCompiledExpression();
    Lepton::CompiledExpression compiled2 = exp2.createCompiledExpression();
    if (compiled1.getVariables().find("x") != compiled1.getVariables().end())
        compiled1.getVariableReference("x") = x;
    if (compiled1.getVariables().find("y") != compiled1.getVariables().end())
        compiled1.getVariableReference("y") = y;
    if (compiled2.getVariables().find("x") != compiled2.getVariables().end())
        compiled2.getVariableReference("x") = x;
    if (compiled2.getVariables().find("y") != compiled2.getVariables().end())
        compiled2.getVariableReference("y") = y;
    assertNumbersEqual(val1, compiled1.evaluate());
    assertNumbersEqual(val2, compiled2.evaluate());

    // Now check CompiledVectorizedExpressions.

    for (int width : Lepton::CompiledVectorExpression::getAllowedWidths()) {
        Lepton::CompiledVectorExpression vector1 = exp1.createCompiledVectorExpression(width);
        Lepton::CompiledVectorExpression vector2 = exp2.createCompiledVectorExpression(width);
        for (int i = 0; i < width; i++) {
            if (vector1.getVariables().find("x") != vector1.getVariables().end())
                for (int j = 0; j < width; j++)
                    vector1.getVariablePointer("x")[j] = (i == j ? x : -100.0);
            if (vector1.getVariables().find("y") != vector1.getVariables().end())
                for (int j = 0; j < width; j++)
                    vector1.getVariablePointer("y")[j] = (i == j ? y : -100.0);
            if (vector2.getVariables().find("x") != vector2.getVariables().end())
                for (int j = 0; j < width; j++)
                    vector2.getVariablePointer("x")[j] = (i == j ? x : -100.0);
            if (vector2.getVariables().find("y") != vector2.getVariables().end())
                for (int j = 0; j < width; j++)
                    vector2.getVariablePointer("y")[j] = (i == j ? y : -100.0);
            const float *result1 = vector1.evaluate();
            const float *result2 = vector2.evaluate();
            assertNumbersEqual(val1, result1[i], 5e-6);
            assertNumbersEqual(val2, result2[i], 5e-6);
        }
    }
}

/**
 * Verify that the derivative of an expression is calculated correctly.
 */

void verifyDerivative(const std::string &expression, const std::string &expectedDeriv)
{
    if (verbose) std::cout << "Checking derivative of: " << expression << "\n";
    Lepton::ParsedExpression computed =
        Lepton::Parser::parse(expression).differentiate("x").optimize();
    Lepton::ParsedExpression expected = Lepton::Parser::parse(expectedDeriv);
    verifySameValue(computed, expected, 1.0, 2.0);
    verifySameValue(computed, expected, 2.0, 3.0);
    verifySameValue(computed, expected, -2.0, 3.0);
    verifySameValue(computed, expected, 2.0, -3.0);
    verifySameValue(computed, expected, 0.0, -3.0);
    verifySameValue(computed, expected, 2.0, 0.0);
}

/**
 * Test the use of a custom function.
 */

void testCustomFunction(const std::string &expression, const std::string &equivalent)
{
    if (verbose) std::cout << "Checking custom function expression: " << expression << "\n";
    std::map<std::string, Lepton::CustomFunction *> functions;
    ExampleFunction exp;
    functions["custom"]           = &exp;
    Lepton::ParsedExpression exp1 = Lepton::Parser::parse(expression, functions);
    Lepton::ParsedExpression exp2 = Lepton::Parser::parse(equivalent);
    verifySameValue(exp1, exp2, 1.0, 2.0);
    verifySameValue(exp1, exp2, 2.0, 3.0);
    verifySameValue(exp1, exp2, -2.0, 3.0);
    verifySameValue(exp1, exp2, 2.0, -3.0);
    Lepton::ParsedExpression deriv1 = exp1.differentiate("x").optimize();
    Lepton::ParsedExpression deriv2 = exp2.differentiate("x").optimize();
    verifySameValue(deriv1, deriv2, 1.0, 2.0);
    verifySameValue(deriv1, deriv2, 2.0, 3.0);
    verifySameValue(deriv1, deriv2, -2.0, 3.0);
    verifySameValue(deriv1, deriv2, 2.0, -3.0);
    Lepton::ParsedExpression deriv3 = deriv1.differentiate("y").optimize();
    Lepton::ParsedExpression deriv4 = deriv2.differentiate("y").optimize();
    verifySameValue(deriv3, deriv4, 1.0, 2.0);
    verifySameValue(deriv3, deriv4, 2.0, 3.0);
    verifySameValue(deriv3, deriv4, -2.0, 3.0);
    verifySameValue(deriv3, deriv4, 2.0, -3.0);
}

TEST(Lepton, Evaluation)
{
    verifyEvaluation("5", 5.0);
    verifyEvaluation("5*2", 10.0);
    verifyEvaluation("2*3+4*5", 26.0);
    verifyEvaluation("2^-3", 0.125);
    verifyEvaluation("1e+2", 100.0);
    verifyEvaluation("-x", 2.0, 3.0, -2.0);
    verifyEvaluation("y^-x", 3.0, 2.0, 0.125);
    verifyEvaluation("1/-x", 3.0, 2.0, -1.0 / 3.0);
    verifyEvaluation("2.1e-4*x*(y+1)", 3.0, 1.0, 1.26e-3);
    verifyEvaluation("sin(2.5)", std::sin(2.5));
    verifyEvaluation("cot(x)", 3.0, 1.0, 1.0 / std::tan(3.0));
    verifyEvaluation("log(x)", 3.0, 1.0, std::log(3.0));
    verifyEvaluation("x^2+y^3+x^-1+y^(1/2)", 1.0, 1.0, 4.0);
    verifyEvaluation("(2*x)*3", 4.0, 4.0, 24.0);
    verifyEvaluation("(x*2)*3", 4.0, 4.0, 24.0);
    verifyEvaluation("2*(x*3)", 4.0, 4.0, 24.0);
    verifyEvaluation("2*(3*x)", 4.0, 4.0, 24.0);
    verifyEvaluation("2*x/3", 1.0, 4.0, 2.0 / 3.0);
    verifyEvaluation("x*2/3", 1.0, 4.0, 2.0 / 3.0);
    verifyEvaluation("5*(-x)*(-y)", 1.0, 4.0, 20.0);
    verifyEvaluation("5*(-x)*(y)", 1.0, 4.0, -20.0);
    verifyEvaluation("5*(x)*(-y)", 1.0, 4.0, -20.0);
    verifyEvaluation("5*(-x)/(-y)", 1.0, 4.0, 1.25);
    verifyEvaluation("5*(-x)/(y)", 1.0, 4.0, -1.25);
    verifyEvaluation("5*(x)/(-y)", 1.0, 4.0, -1.25);
    verifyEvaluation("x+(-y)", 1.0, 4.0, -3.0);
    verifyEvaluation("(-x)+y", 1.0, 4.0, 3.0);
    verifyEvaluation("x/(1/y)", 1.0, 4.0, 4.0);
    verifyEvaluation("x*w; w = 5", 3.0, 1.0, 15.0);
    verifyEvaluation("a+b^2;a=x-b;b=3*y", 2.0, 3.0, 74.0);
    verifyEvaluation("erf(x)+erfc(x)", 2.0, 3.0, 1.0);
    verifyEvaluation("min(3, x)", 2.0, 3.0, 2.0);
    verifyEvaluation("min(y, 5)", 2.0, 3.0, 3.0);
    verifyEvaluation("max(x, y)", 2.0, 3.0, 3.0);
    verifyEvaluation("max(x, -1)", 2.0, 3.0, 2.0);
    verifyEvaluation("abs(x-y)", 2.0, 3.0, 1.0);
    verifyEvaluation("delta(x)+3*delta(y-1.5)", 2.0, 1.5, 3.0);
    verifyEvaluation("step(x-3)+y*step(x)", 2.0, 3.0, 3.0);
    verifyEvaluation("floor(x)", -2.1, 3.0, -3.0);
    verifyEvaluation("ceil(x)", -2.1, 3.0, -2.0);
    verifyEvaluation("select(x, 1.0, y)", 0.3, 2.0, 1.0);
    verifyEvaluation("select(x, 1.0, y)", 0.0, 2.0, 2.0);
    verifyEvaluation("atan2(x, y)", 3.0, 1.5, std::atan(2.0));
    verifyEvaluation("sqrt(x^2)", -2.2, 0.0, 2.2);
    verifyEvaluation("sqrt(x)^2", 2.2, 0.0, 2.2);
    verifyEvaluation("x^2+x^4", 2.0, 0.0, 20.0);
    verifyEvaluation("x^-2+x^-3", 2.0, 0.0, 0.375);
    verifyEvaluation("x^1.8", 2.2, 0.0, std::pow(2.2, 1.8));
}

TEST(Lepton, InvalidEvaluation)
{
    ASSERT_NO_THROW(verifyInvalidExpression("1..2"));
    ASSERT_NO_THROW(verifyInvalidExpression("1*(2+3"));
    ASSERT_NO_THROW(verifyInvalidExpression("5++4"));
    ASSERT_NO_THROW(verifyInvalidExpression("1+2)"));
    ASSERT_NO_THROW(verifyInvalidExpression("cos(2,3)"));
}

TEST(Lepton, VerifyDerivative)
{
    verifyDerivative("x", "1");
    verifyDerivative("x^2+x", "2*x+1");
    verifyDerivative("y^x-x", "log(y)*(y^x)-1");
    verifyDerivative("sin(x)", "cos(x)");
    verifyDerivative("cos(x)", "-sin(x)");
    verifyDerivative("tan(x)", "square(sec(x))");
    verifyDerivative("cot(x)", "-square(csc(x))");
    verifyDerivative("sec(x)", "sec(x)*tan(x)");
    verifyDerivative("csc(x)", "-csc(x)*cot(x)");
    verifyDerivative("exp(2*x)", "2*exp(2*x)");
    verifyDerivative("log(x)", "1/x");
    verifyDerivative("sqrt(x)", "0.5/sqrt(x)");
    verifyDerivative("asin(x)", "1/sqrt(1-x^2)");
    verifyDerivative("acos(x)", "-1/sqrt(1-x^2)");
    verifyDerivative("atan(x)", "1/(1+x^2)");
    verifyDerivative("atan2(2*x,y)", "2*y/(4*x^2+y^2)");
    verifyDerivative("sinh(x)", "cosh(x)");
    verifyDerivative("cosh(x)", "sinh(x)");
    verifyDerivative("tanh(x)", "1/(cosh(x)^2)");
    verifyDerivative("erf(x)", "1.12837916709551*exp(-x^2)");
    verifyDerivative("erfc(x)", "-1.12837916709551*exp(-x^2)");
    verifyDerivative("step(x)*x+step(1-x)*2*x", "step(x)+step(1-x)*2");
    verifyDerivative("recip(x)", "-1/x^2");
    verifyDerivative("square(x)", "2*x");
    verifyDerivative("cube(x)", "3*x^2");
    verifyDerivative("min(x, 2*x)", "step(x-2*x)*2+(1-step(x-2*x))*1");
    verifyDerivative("max(5, x^2)", "(1-step(5-x^2))*2*x");
    verifyDerivative("abs(3*x)", "step(3*x)*3+(1-step(3*x))*-3");
    verifyDerivative("floor(x)+0.5*x*ceil(x)", "0.5*ceil(x)");
    verifyDerivative("select(x, x^2, 3*x)", "select(x, 2*x, 3)");
}

TEST(Lepton, CustomFunction)
{
    testCustomFunction("custom(x, y)/2", "x*y");
    testCustomFunction("custom(x^2, 1)+custom(2, y-1)", "2*x^2+4*(y-1)");
}

TEST(Lepton, Optimize)
{
    std::string buffer;
    std::stringstream out(buffer);

    out << Lepton::Parser::parse("x*x").optimize();
    ASSERT_THAT(out.str(), StrEq("square(x)"));
    out.str("");

    out << Lepton::Parser::parse("x*x*x").optimize();
    ASSERT_THAT(out.str(), StrEq("cube(x)"));
    out.str("");

    out << Lepton::Parser::parse("x*(x*x)").optimize();
    ASSERT_THAT(out.str(), StrEq("cube(x)"));
    out.str("");

    out << Lepton::Parser::parse("(x*x)*x").optimize();
    ASSERT_THAT(out.str(), StrEq("cube(x)"));
    out.str("");

    out << Lepton::Parser::parse("2*3*x").optimize();
    ASSERT_THAT(out.str(), StrEq("6*(x)"));
    out.str("");

    out << Lepton::Parser::parse("1/(1+x)").optimize();
    ASSERT_THAT(out.str(), StrEq("recip(1+(x))"));
    out.str("");

    out << Lepton::Parser::parse("x^(1/2)").optimize();
    ASSERT_THAT(out.str(), StrEq("sqrt(x)"));
    out.str("");
    out << Lepton::Parser::parse("log(3*cos(x))^(sqrt(4)-2)").optimize();
    ASSERT_THAT(out.str(), StrEq("1"));
    out.str("");
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }
    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
