// nfsim_funcparser.h — ExprTk-based drop-in replacement for mu::Parser in NFsim
//
// Maintained by bngsim to keep expression handling aligned between the
// ODE/SSA engines and the vendored NFsim runtime.
//
// This header provides a `mu::Parser` class with the same API surface as muParser
// but backed by ExprTk. All existing NFsim code that uses `mu::Parser` works without
// changes to variable names, namespace references, or calling patterns.
//
// Key design decisions:
// - "Constants" are stored as mutable variables in internal storage, so
//   DefineConst() can be called after SetExpr() (used by updateParameters/fileUpdate).
// - ExprTk's `log` is natural log (matching bngsim). BNG's `ln` is aliased to natural log.
// - Built-in functions match bngsim's expression.cpp: if(), ln(), rint(), sign(), etc.
// - Constants: _PI, _e, _Na (NFsim convention) plus bngsim's _pi, _kB, _NA, _R, _h, _F.
//
// IMPORTANT: ExprTk does not allow identifiers starting with underscore ('_').
// NFsim uses _PI, _e, _Na as constant names, and BNG XML files may define
// parameters with leading underscores (e.g., __FREE parameters in PyBNF,
// __TFUN_VAL__ placeholders for time-dependent functions).
// This shim transparently remaps: "_X" → "u_X" in both symbol registration
// and expression preprocessing. The caller never sees this — they continue to
// use DefineConst("_PI", ...) and SetExpr("sin(_e * _PI)") as before.
//
#ifndef NFSIM_FUNCPARSER_H_
#define NFSIM_FUNCPARSER_H_

// ExprTk compilation options — disable features we don't need.
#define exprtk_disable_string_capabilities
#define exprtk_disable_rtl_io_file
#define exprtk_disable_rtl_vecops
// BNG is case-sensitive for parameter names (e.g., k3 ≠ K3).
// ExprTk defaults to case-insensitive, which silently merges k3/K3.
#define exprtk_disable_caseinsensitivity
#include "exprtk.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace mu {

// ─── Custom ExprTk function adapters ─────────────────────────────────────────

namespace detail {

// 1-arg: ln(x) — natural logarithm (backward-compat alias)
template <typename T>
struct LnFunction : public exprtk::ifunction<T> {
    LnFunction() : exprtk::ifunction<T>(1) {}
    T operator()(const T& x) override { return std::log(x); }
};

// 1-arg: rint(x) — round to nearest integer as muParser's rint does,
// floor(x + 0.5): a half rounds up, so rint(-2.5) is -2 and rint(-0.5) is 0.
// BioNetGen's own rint is the same (BNG2.pl, run_network). std::round rounds a
// half away from zero and so differed at every negative half.
template <typename T>
struct RintFunction : public exprtk::ifunction<T> {
    RintFunction() : exprtk::ifunction<T>(1) {}
    T operator()(const T& x) override { return std::floor(x + T(0.5)); }
};

// 1-arg: sign(x) — signum function
template <typename T>
struct SignFunction : public exprtk::ifunction<T> {
    SignFunction() : exprtk::ifunction<T>(1) {}
    T operator()(const T& x) override {
        return (x > T(0)) ? T(1) : ((x < T(0)) ? T(-1) : T(0));
    }
};

// 0-arg: a function the caller defines with DefineFun(name, f).
template <typename T>
struct NullaryFunction : public exprtk::ifunction<T> {
    explicit NullaryFunction(T (*f)()) : exprtk::ifunction<T>(0), f_(f) {}
    T operator()() override { return f_(); }
    T (*f_)();
};

// ─── Underscore name remapping ───────────────────────────────────────────────
// ExprTk rejects identifiers starting with '_'. NFsim/BNG uses _PI, _e, _Na,
// and TFUN injects __TFUN_VAL__ after SetExpr(). We remap: "_X" → "u_X"
// transparently in both symbol names and expressions.

inline std::string remap_name(const std::string& name) {
    if (!name.empty() && name[0] == '_') {
        return "u_" + name.substr(1);
    }
    return name;
}

// Remap ALL underscore-prefixed identifiers in an expression string.
// This scans the expression for any token starting with '_' and applies
// remap_name() to it, regardless of whether it was previously registered
// via DefineConst/DefineVar. This is critical for TFUN support where
// __TFUN_VAL__ is injected via DefineConst *after* SetExpr() — the
// expression string must be remapped at SetExpr() time even though the
// name hasn't been tracked yet.
inline std::string remap_expression(const std::string& expr,
    const std::vector<std::string>& /* underscore_names — kept for API compat */)
{
    std::string result;
    result.reserve(expr.size() + 16);
    size_t i = 0;
    while (i < expr.size()) {
        // Check for underscore-leading identifier at a word boundary
        if (expr[i] == '_') {
            // Verify it's at a word boundary (not mid-identifier)
            bool at_boundary = (i == 0) ||
                (!std::isalnum(expr[i - 1]) && expr[i - 1] != '_');
            if (at_boundary) {
                // Collect the full identifier: _[A-Za-z0-9_]+
                size_t start = i;
                i++;  // skip the leading '_'
                while (i < expr.size() &&
                       (std::isalnum(expr[i]) || expr[i] == '_')) {
                    i++;
                }
                std::string token = expr.substr(start, i - start);
                result += remap_name(token);
                continue;
            }
        }
        result += expr[i];
        i++;
    }
    return result;
}

// ─── Numeric literals, read as muParser reads them ──────────────────────────
// muParser reads a number with a C++ stream (stringstream >> double), which
// gives the double nearest the decimal, as strtod does. ExprTk's own
// string_to_real sums the digits in a double and then scales by powers of ten,
// rounding at each step, so it is not correctly rounded: 7.77e-05 came back
// 7.769999999999999e-05, one unit in the last place low, and
// 0.49999999999999997 came back 0.5. exact_literals() respells each literal
// ExprTk would misread as a form it reads exactly. A literal ExprTk reads
// right is left as written, and so is everything else.

// End of the number token ExprTk's lexer scans from i (scan_number), or npos
// where the lexer would report a malformed number.
inline std::size_t number_token_end(const std::string& expr, std::size_t i) {
    bool dot = false, exp = false, exp_sign = false, exp_digit = false;
    const auto digit = [](char c) { return std::isdigit(static_cast<unsigned char>(c)) != 0; };
    while (i < expr.size()) {
        const char c = expr[i];
        if (c == '.') {
            if (dot) return std::string::npos;
            dot = true;
        } else if (c == 'e' || c == 'E') {
            if (i + 1 >= expr.size() ||
                (expr[i + 1] != '+' && expr[i + 1] != '-' && !digit(expr[i + 1]))) {
                return std::string::npos;
            }
            exp = true;
        } else if (exp && (c == '+' || c == '-') && !exp_digit) {
            if (exp_sign) return std::string::npos;
            exp_sign = true;
        } else if (exp && digit(c)) {
            exp_digit = true;
        } else if (!digit(c)) {
            break;
        }
        ++i;
    }
    return i;
}

// End of the symbol token ExprTk's lexer scans from i, a letter or an
// underscore (scan_symbol): word characters, and a '.' another word character
// follows.
inline std::size_t symbol_token_end(const std::string& expr, std::size_t i) {
    const auto word = [](char c) {
        return std::isalnum(static_cast<unsigned char>(c)) != 0 || c == '_';
    };
    std::size_t j = i + 1;
    while (j < expr.size() &&
           (word(expr[j]) || (expr[j] == '.' && j + 1 < expr.size() && word(expr[j + 1])))) {
        ++j;
    }
    return j;
}

inline bool exprtk_reads(const std::string& text, double& value) {
    return exprtk::details::string_to_real(text, value);
}

// A spelling of the finite, non-negative double v that ExprTk reads as v.
inline std::string exact_spelling(double v) {
    // v's shortest round-trip digits as an integer D and an exponent E. ExprTk
    // reads an integer below 2^53 exactly and scales it by one exactly
    // representable power of ten in one multiply or divide, which rounds once,
    // to the double nearest D*10^E: strtod's answer, and v. Checked rather
    // than assumed, since a 16- or 17-digit D or a larger exponent falls
    // outside that.
    for (int digits = 1; digits <= 17; ++digits) {
        char buf[40];
        std::snprintf(buf, sizeof(buf), "%.*e", digits - 1, v);
        if (std::strtod(buf, nullptr) != v) continue;
        const char* e = std::strchr(buf, 'e');
        std::string mantissa;
        for (const char* c = buf; c != e; ++c) {
            if (*c != '.') mantissa += *c;
        }
        int exponent = std::atoi(e + 1) - (digits - 1);
        while (mantissa.size() > 1 && mantissa.back() == '0') {
            mantissa.pop_back();
            ++exponent;
        }
        const std::string candidate =
            exponent == 0 ? mantissa : mantissa + "e" + std::to_string(exponent);
        double read = 0.0;
        if (exprtk_reads(candidate, read) && read == v) return candidate;
        break;
    }
    // Otherwise exactly: v is M*2^k with M an odd integer below 2^53. ExprTk
    // reads each factor exactly, and a multiply or divide by a power of two is
    // exact, so the constant it folds is v itself.
    int exp2 = 0;
    const double frac = std::frexp(v, &exp2);
    long long mant = static_cast<long long>(std::ldexp(frac, 53));
    int scale = exp2 - 53;
    while (mant != 0 && mant % 2 == 0) {
        mant /= 2;
        ++scale;
    }
    std::string out = "(" + std::to_string(mant);
    for (int left = std::abs(scale); left > 0;) {
        const int step = std::min(left, 52);
        out += (scale < 0 ? '/' : '*') + std::to_string(1LL << step);
        left -= step;
    }
    return out + ")";
}

inline std::string exact_literals(const std::string& expr) {
    const auto letter = [](char c) { return std::isalpha(static_cast<unsigned char>(c)) != 0; };
    const auto digit = [](char c) { return std::isdigit(static_cast<unsigned char>(c)) != 0; };
    std::string out;
    out.reserve(expr.size());
    std::size_t i = 0;
    while (i < expr.size()) {
        const char c = expr[i];
        if (letter(c) || c == '_') {
            // A symbol, copied whole so that its digits are never read as a number.
            const std::size_t j = symbol_token_end(expr, i);
            out.append(expr, i, j - i);
            i = j;
            continue;
        }
        if (digit(c) || (c == '.' && i + 1 < expr.size() && digit(expr[i + 1]))) {
            const std::size_t j = number_token_end(expr, i);
            if (j == std::string::npos) {
                out.append(expr, i, std::string::npos);  // ExprTk reports it
                break;
            }
            const std::string token = expr.substr(i, j - i);
            char* end = nullptr;
            const double v = std::strtod(token.c_str(), &end);
            double read = 0.0;
            // A literal ExprTk cannot read at all (1e-320, an exponent below
            // its range) is left for ExprTk to report.
            if (end == token.c_str() + token.size() && std::isfinite(v) &&
                exprtk_reads(token, read) && read != v) {
                const std::string spelled = exact_spelling(v);
                // ExprTk multiplies a closing bracket by a number after it,
                // but a closing bracket followed by a bracket is a second
                // statement, so the product is written out.
                const std::size_t last = out.find_last_not_of(" \t\n\v\f\r");
                if (spelled[0] == '(' && last != std::string::npos &&
                    (out[last] == ')' || out[last] == ']' || out[last] == '}')) {
                    out += '*';
                }
                out += spelled;
            } else {
                out += token;
            }
            i = j;
            continue;
        }
        out += c;
        ++i;
    }
    return out;
}

}  // namespace detail

// ─── mu::Parser — ExprTk-based drop-in replacement ──────────────────────────

class Parser {
public:
    struct exception_type : public std::runtime_error {
        exception_type(const std::string& msg) : std::runtime_error(msg) {}
        std::string GetMsg() const { return what(); }
    };

    Parser()
        : compiled_(false)
    {
        // Register built-in constants (NFsim convention: _PI, _e, _Na)
        DefineConst("_PI", 3.14159265358979323846);
        DefineConst("_e",  2.71828182845904523536);
        DefineConst("_Na", 6.02214076e23);

        // Additional constants supported by bngsim: _pi, _NA, _kB, _R, _h, _F
        DefineConst("_pi", 3.14159265358979323846);
        DefineConst("_NA", 6.02214076e23);
        DefineConst("_kB", 1.380649e-23);
        DefineConst("_R",  8.314462618153241);
        DefineConst("_h",  6.62607015e-34);
        DefineConst("_F",  96485.33212331002);

        // Register backward-compatible aliases
        symbol_table_.add_function("ln",   ln_func_);
        symbol_table_.add_function("rint", rint_func_);
        symbol_table_.add_function("sign", sign_func_);
    }

    ~Parser() = default;

    // Non-copyable (symbol_table holds pointers to internal storage)
    Parser(const Parser&) = delete;
    Parser& operator=(const Parser&) = delete;

    // ─── DefineVar: bind a variable name to an external double* ─────────
    void DefineVar(const std::string& name, double* ptr) {
        std::string mapped = detail::remap_name(name);
        if (name != mapped) {
            // Track underscore-prefixed name for expression remapping
            trackUnderscoreName(name);
        }
        if (compiled_) {
            symbol_table_.add_variable(mapped, *ptr);
            recompile();
        } else {
            symbol_table_.add_variable(mapped, *ptr);
        }
    }

    // ─── DefineConst: store value internally, register as variable ──────
    void DefineConst(const std::string& name, double value) {
        std::string mapped = detail::remap_name(name);
        if (name != mapped) {
            trackUnderscoreName(name);
        }
        auto it = const_storage_.find(mapped);
        if (it != const_storage_.end()) {
            *(it->second) = value;
        } else {
            auto ptr = std::make_unique<double>(value);
            double* raw = ptr.get();
            const_storage_[mapped] = std::move(ptr);
            symbol_table_.add_variable(mapped, *raw);
        }
    }

    // ─── DefineFun: a function of no arguments, as muParser defines one ──
    void DefineFun(const std::string& name, double (*f)()) {
        nullary_funcs_.push_back(std::make_unique<detail::NullaryFunction<double>>(f));
        symbol_table_.add_function(detail::remap_name(name), *nullary_funcs_.back());
    }

    // ─── SetExpr: compile the expression ────────────────────────────────
    void SetExpr(const std::string& expr) {
        original_expr_string_ = expr;
        // Remap underscore-prefixed identifiers before ExprTk compilation
        expr_string_ = detail::remap_expression(expr, underscore_names_);
        // Respell any literal ExprTk would read otherwise than muParser does.
        // That is to change a value, never whether the expression is accepted,
        // so the expression as written is compiled first: one ExprTk refuses
        // (two numbers side by side, say) stays refused, where a literal
        // respelled in brackets would join its neighbour as a product.
        const std::string exact = detail::exact_literals(expr_string_);
        if (exact != expr_string_) {
            compile();
            expr_string_ = exact;
        }
        compile();
    }

    // ─── Eval: evaluate the compiled expression ─────────────────────────
    double Eval() {
        if (!compiled_) {
            throw exception_type("Expression not compiled (call SetExpr first)");
        }
        return expression_.value();
    }

    // ─── GetExpr: return the ORIGINAL expression string ─────────────────
    std::string GetExpr() const {
        return original_expr_string_;
    }

private:
    void trackUnderscoreName(const std::string& name) {
        for (const auto& n : underscore_names_) {
            if (n == name) return;  // already tracked
        }
        underscore_names_.push_back(name);
    }

    void compile() {
        expression_ = exprtk::expression<double>();  // reset
        expression_.register_symbol_table(symbol_table_);

        exprtk::parser<double> parser;
        // Increase max stack depth for deeply nested if() expressions.
        // ExprTk default is 400 (~200 nested if()), muParser handled 2000.
        parser.settings().set_max_stack_depth(4096);
        if (!parser.compile(expr_string_, expression_)) {
            throw exception_type(
                "ExprTk compilation failed for '" + original_expr_string_ +
                "' (remapped: '" + expr_string_ + "'): " +
                parser.error());
        }
        compiled_ = true;
    }

    void recompile() {
        if (!expr_string_.empty()) {
            compile();
        }
    }

    // ExprTk objects
    exprtk::symbol_table<double> symbol_table_;
    exprtk::expression<double> expression_;
    std::string expr_string_;           // remapped expression (what ExprTk sees)
    std::string original_expr_string_;  // original expression (what caller sees)
    bool compiled_;

    // Internal storage for "constants" (mutable via DefineConst after SetExpr)
    std::unordered_map<std::string, std::unique_ptr<double>> const_storage_;

    // Track names that start with '_' for expression remapping
    std::vector<std::string> underscore_names_;

    // Custom function objects (must outlive symbol_table_)
    detail::LnFunction<double> ln_func_;
    detail::RintFunction<double> rint_func_;
    detail::SignFunction<double> sign_func_;
    std::vector<std::unique_ptr<detail::NullaryFunction<double>>> nullary_funcs_;
};

}  // namespace mu

#endif  // NFSIM_FUNCPARSER_H_
