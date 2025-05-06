#!/bin/bash

# Function to check if a container exists
check_container() {
    local container=$1
    echo "Checking container: $container"
    if docker pull $container 2>&1 | grep -q "Digest: sha256:"; then
        echo "✓ Container $container is available"
        # Kill the pull immediately after checking digest
        pkill -9 -f "docker pull $container"
        return 0
    else
        echo "✗ Container $container is not available"
        return 1
    fi
}

# Array of required containers
containers=(
    "quay.io/biocontainers/fastqc:0.11.9--0"
    "quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0"
    "quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2"
    "quay.io/biocontainers/star:2.7.9a--h9ee0642_0"
    "quay.io/biocontainers/salmon:1.5.2--h84f40af_0"
    "quay.io/biocontainers/subread:2.0.1--hed695b0_0"
    "quay.io/biocontainers/bioconductor-deseq2:1.32.0--r41h399db7b_0"
)

# Counter for failed checks
failed=0

echo "Checking availability of required containers..."
echo "----------------------------------------------"

# Check each container
for container in "${containers[@]}"; do
    if ! check_container "$container"; then
        ((failed++))
    fi
    echo "----------------------------------------------"
done

# Final status
echo "Container check completed."
if [ $failed -eq 0 ]; then
    echo "All containers are available! ✓"
    exit 0
else
    echo "Error: $failed container(s) are not available! ✗"
    exit 1
fi 