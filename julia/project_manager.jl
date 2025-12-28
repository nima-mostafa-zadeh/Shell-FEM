# Start Julia in your project directory

# Activate the project environment
using Pkg
Pkg.activate(".")

# Generate Project.toml
#Pkg.generate("EnucleationFEM")

# Or if already in a directory with code:
# The Project.toml will be created when you activate

# Add your dependencies
Pkg.add("Plots")