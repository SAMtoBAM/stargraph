#!/bin/bash

set -e

echo "Installing Geneml into the stargraph environment..."

"$PREFIX/bin/python" -m pip install --no-cache-dir geneml

echo "Geneml installation completed."
