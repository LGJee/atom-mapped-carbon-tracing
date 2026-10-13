#!/usr/bin/env bash
set -euo pipefail

root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
prefix="${1:-$root/.envs}"
if ! command -v conda >/dev/null 2>&1; then
    echo 'Conda was not found. Install Miniconda/Miniforge and initialize this shell first.' >&2
    exit 1
fi
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
    echo 'These pinned environment specifications target Linux/WSL x86_64.' >&2
    exit 1
fi
command -v git >/dev/null
mkdir -p -- "$prefix"
for name in localmapper rxnmapper enzymemap; do
    env="$prefix/$name"
    if [[ -x "$env/bin/python" ]]; then
        echo "Checking existing environment: $env"
    else
        conda env create --prefix "$env" --file "$root/environments/$name.yml"
    fi
    "$env/bin/python" "$root/scripts/check_environment.py" "$name"
done
echo 'Environment setup and real mapper smoke tests completed.'
echo 'If a custom prefix was used, set environment_dir in your JSON configuration.'
