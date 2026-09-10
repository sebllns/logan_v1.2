# Logan v1.2

## SETUP

Environment: https://docs.tacc.utexas.edu/hpc/vista/
Architecture: ARM
TACC ulimit: 2560000

Download and add the `kmindex` binaries in `PATH`:
https://github.com/tlemane/kmindex/releases/download/v0.6.1/kmindex-v0.6.1-Linux-arm64.tar.gz
https://github.com/tlemane/kmtricks/releases/download/v1.6.0/kmtricks-v1.6.0-Linux-arm64.tar.gz


Install and test `kmhelpers`
```bash
module load gcc cuda
module load python3
git clone https://github.com/sebllns/kmhelpers.git
cd kmhelpers
python3.11 -m venv .env
source .env/bin/activate
pip3 install -e .[dev]
pytest
```

Use kmhelpers
```bash
module load gcc cuda
module load python3
cd kmhelpers
source .env/bin/activate
kmhelpers --version
```
