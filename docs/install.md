# Install

LAMPpy is pure Python and needs **Python 3.10 or newer**.

```bash
pip install lamppy
```

Development install, from a clone:

```bash
pip install -e ".[dev]"
```

## External tools

The Python install gives you `kmers`, `targets`, the scoring, the report and the
thermodynamics. The stages that delegate need their tools on `PATH`:

| Stage | Needs | Why |
|---|---|---|
| `design --backend primer3` | nothing | designs in Python with `primer3-py` |
| `design --backend glapd` | [GLAPD](https://github.com/jiqingxiaoxi/GLAPD) (built with `make`), `perl`, `bowtie`, `bowtie-build` | primer design and the specificity filter |
| `validate` | BLAST+ (`blastn`) | short-query site discovery |
| `prepare --dustmasker` | BLAST+ (`dustmasker`) | low-complexity soft-masking |
| `prepare --trf` | [`trf`](https://tandem.bu.edu/trf/trf.html) | tandem-repeat masking |

None of these is on PyPI, so `pip` cannot install them. Everything except GLAPD is on
bioconda, pinned in `environment.yml`:

```bash
mamba env create -f environment.yml
mamba activate lamppy
pip install -e ".[dev]"
```

GLAPD has no package anywhere. Build it once and point `GLAPD` at the directory, which
`design --glapd` then reads by default:

```bash
git clone https://github.com/jiqingxiaoxi/GLAPD.git ~/GLAPD && make -C ~/GLAPD
export GLAPD=~/GLAPD
```

Check what is present, GLAPD included:

```bash
lamppy check-deps
```

## Container

The image pins every external tool and builds GLAPD from a pinned commit, which is the
quickest way to get a working `design` and `validate` on macOS:

```bash
docker build -t lamppy .
docker run --rm -v "$PWD:/data" lamppy check-deps
docker run --rm -v "$PWD:/data" lamppy run --config /data/run.yaml
```

The GLAPD commit the image was built from is recorded in `/opt/GLAPD/COMMIT`.

GLAPD is a Linux-oriented GPL-2.0 source tarball with no conda recipe. LAMPpy invokes its
binaries as separate processes and never vendors or patches its source, which keeps the
licences separate. The bundled `bowtie` binaries in the GLAPD tarball are Linux ELF files;
on macOS install `bowtie` separately (for example `brew install bowtie`).
