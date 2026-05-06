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

#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace mu {

// ─── Custom ExprTk function adapters ─────────────────────────────────────────

namespace detail {

// 3-arg: if(cond, true_val, false_val)
// NOTE: This custom IfFunction is effectively dead code. ExprTk has a
// built-in `if` keyword that takes precedence over add_function("if", ...).
// The built-in uses nonzero truthiness (cond != 0), matching BNG Perl semantics.
template <typename T>
struct IfFunction : public exprtk::ifunction<T> {
    IfFunction() : exprtk::ifunction<T>(3) {
        exprtk::ifunction<T>::allow_zero_parameters() = false;
    }
    T operator()(const T& cond, const T& true_val, const T& false_val) override {
        return (cond > T(0.5)) ? true_val : false_val;  // DEAD CODE — see note above
    }
};

// 1-arg: ln(x) — natural logarithm (backward-compat alias)
template <typename T>
struct LnFunction : public exprtk::ifunction<T> {
    LnFunction() : exprtk::ifunction<T>(1) {}
    T operator()(const T& x) override { return std::log(x); }
};

// 1-arg: rint(x) — round to nearest integer
template <typename T>
struct RintFunction : public exprtk::ifunction<T> {
    RintFunction() : exprtk::ifunction<T>(1) {}
    T operator()(const T& x) override { return std::round(x); }
};

// 1-arg: sign(x) — signum function
template <typename T>
struct SignFunction : public exprtk::ifunction<T> {
    SignFunction() : exprtk::ifunction<T>(1) {}
    T operator()(const T& x) override {
        return (x > T(0)) ? T(1) : ((x < T(0)) ? T(-1) : T(0));
    }
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
        symbol_table_.add_function("if",   if_func_);
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

    // ─── SetExpr: compile the expression ────────────────────────────────
    void SetExpr(const std::string& expr) {
        original_expr_string_ = expr;
        // Remap underscore-prefixed identifiers before ExprTk compilation
        expr_string_ = detail::remap_expression(expr, underscore_names_);
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
    detail::IfFunction<double> if_func_;
    detail::LnFunction<double> ln_func_;
    detail::RintFunction<double> rint_func_;
    detail::SignFunction<double> sign_func_;
};

}  // namespace mu

#endif  // NFSIM_FUNCPARSER_H_
