#!/bin/bash

# Exit on error
set -e

# Print commands
set -x

# Change to project root
cd ..

# Check if results directory exists
if [ -d "results" ]; then
    echo "Deleting results directory..."
    rm -rf results
    echo "Results directory deleted successfully!"
else
    echo "No results directory found."
fi

# Return to original directory
cd - > /dev/null 