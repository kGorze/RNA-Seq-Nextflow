#!/bin/bash

set -e  # Exit on error
set -u  # Exit on undefined variable

echo "Running final pre-commit checks..."

# Unset Python environment variables
unset PYTHONHOME
unset PYTHONPATH

# Check if Python 3 is installed
if ! command -v python3 &> /dev/null; then
    echo "ERROR: Python 3 is not installed!"
    exit 1
fi

# Check if required Python tools are installed
check_python_tool() {
    if ! python3 -m $1 --help &> /dev/null; then
        echo "ERROR: $1 is not installed!"
        echo "Please install it using: pip3 install $1"
        return 1
    fi
    return 0
}

# Check for required Python tools
echo "Checking for required Python tools..."
required_tools=("flake8" "black")
tools_missing=0

for tool in "${required_tools[@]}"; do
    if ! check_python_tool "$tool"; then
        tools_missing=1
    fi
done

if [ $tools_missing -eq 1 ]; then
    echo "Please install missing tools and try again."
    echo "You can install all required tools with:"
    echo "pip3 install flake8 black"
    exit 1
fi

# Check Python code style
echo "Checking Python code style..."
cd .. && {
    python3 -m flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics
    python3 -m black --check .
}
cd - > /dev/null

# Check Nextflow syntax
echo "Checking Nextflow syntax..."
cd .. && {
    nextflow run main.nf --help > /dev/null || {
        echo "ERROR: Nextflow syntax check failed!"
        exit 1
    }
    echo "✓ Nextflow syntax check passed"
}
cd - > /dev/null

# Run minimal test
echo "Running minimal pipeline test..."
./test_pipeline_mini.sh

# Check Docker containers
echo "Checking Docker containers..."
./check_containers.sh
if [ $? -ne 0 ]; then
    echo "ERROR: Container check failed!"
    exit 1
fi

# Check for required files
echo "Checking for required files..."
required_files=(
    "main.nf"
    "nextflow.config"
    "README.md"
    "test_data/test_data_generator.py"
)

for file in "${required_files[@]}"; do
    if [ ! -f "../$file" ]; then
        echo "ERROR: Required file $file not found!"
        exit 1
    fi
done

echo "All checks passed successfully!"
