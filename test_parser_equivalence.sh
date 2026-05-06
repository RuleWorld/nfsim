#!/bin/bash
# Test script to verify muParser and ExprTk backends produce equivalent results
#
# This script builds NFsim with both parsers and runs a simple model to ensure
# they produce identical output.

set -e

echo "=== Testing Parser Equivalence ==="
echo ""

# Clean any existing builds
rm -rf build_muparser build_exprtk

# Build with muParser (default)
echo "Building with muParser..."
mkdir -p build_muparser
cd build_muparser
cmake -DNFSIM_USE_EXPRTK=OFF .. > /dev/null
make -j4 > /dev/null 2>&1
cd ..
echo "  ✓ muParser build successful"

# Build with ExprTk
echo "Building with ExprTk..."
mkdir -p build_exprtk
cd build_exprtk
cmake -DNFSIM_USE_EXPRTK=ON .. > /dev/null
make -j4 > /dev/null 2>&1
cd ..
echo "  ✓ ExprTk build successful"

# Find a simple test model
TEST_MODEL="models/simple_system/simple_system.xml"
if [ ! -f "$TEST_MODEL" ]; then
    echo "Error: Test model not found at $TEST_MODEL"
    exit 1
fi

# Run with both parsers
echo ""
echo "Running test model with both parsers..."
build_muparser/NFsim -xml "$TEST_MODEL" -o muparser_output.gdat -seed 12345 -sim 10 > /dev/null 2>&1
build_exprtk/NFsim -xml "$TEST_MODEL" -o exprtk_output.gdat -seed 12345 -sim 10 > /dev/null 2>&1

# Compare outputs
if diff -q muparser_output.gdat exprtk_output.gdat > /dev/null 2>&1; then
    echo "  ✓ Outputs match! Parsers are equivalent."
    rm muparser_output.gdat exprtk_output.gdat
    echo ""
    echo "=== All tests passed ==="
    exit 0
else
    echo "  ✗ Outputs differ! Parsers produce different results."
    echo ""
    echo "Diff:"
    diff muparser_output.gdat exprtk_output.gdat || true
    exit 1
fi
