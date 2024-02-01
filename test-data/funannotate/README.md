# funannotate test-data

## funnannotate databases

You need to specify AUGUSTUS_CONFIG_PATH to install the databases.
It's present at `/usr/local/config/` inside the container.1
Seems to work, but takes ages.

```bash
apptainer exec \
    -B ${PWD},${TMPDIR} \
    -H $(mktemp -d) \
    --pwd ${PWD} \
    --containall \
    --cleanenv \
    --writable-tmpfs \
    --env AUGUSTUS_CONFIG_PATH=/usr/local/config \
    docker://ghcr.io/tomharrop/container-funannotate:1.8.15_cv2 \
        funannotate setup \
            -i all \
            -d test-data/funannotate/db
```

## get the augustus DB from the container

It needs to be writable (!) to run some of the tools.

```bash
tmpdir="$( mktemp -d )"
AUGUSTUS_CONFIG_PATH="$( readlink -f test-data/funannotate/augustus_config )"
mkdir -p "$( dirname "${AUGUSTUS_CONFIG_PATH}" )"

apptainer exec \
    -B ${tmpdir} \
    -H ${tmpdir} \
    --pwd ${tmpdir} \
    --containall \
    --cleanenv \
    --writable-tmpfs \
    docker://ghcr.io/tomharrop/container-funannotate:1.8.15_cv3 \
        cp -r /usr/local/config ./

mv "${tmpdir}/config" "${AUGUSTUS_CONFIG_PATH}"

```

Then run funannotate like this:

```bash
apptainer exec \
    -B ${PWD},${TMPDIR},${AUGUSTUS_CONFIG_PATH} \
    -H $(mktemp -d) \
    --pwd ${PWD} \
    --containall \
    --cleanenv \
    --writable-tmpfs \
    --env AUGUSTUS_CONFIG_PATH="${AUGUSTUS_CONFIG_PATH}" \
    docker://ghcr.io/tomharrop/container-funannotate:1.8.15_cv3 \
        funannotate setup \
            -i all \
            -d test-data/funannotate/db

apptainer exec \
    -B ${PWD},${TMPDIR},${AUGUSTUS_CONFIG_PATH} \
    -H $(mktemp -d) \
    --pwd ${PWD} \
    --containall \
    --cleanenv \
    --writable-tmpfs \
    --env AUGUSTUS_CONFIG_PATH="${AUGUSTUS_CONFIG_PATH}" \
    --env FUNANNOTATE_DB=test-data/funannotate/db \
    docker://ghcr.io/tomharrop/container-funannotate:1.8.15_cv3 \
        funannotate setup \
            -b embryophyta
```


embryophyta

## get the eggnog DB

```bash

mkdir -p test-data/funannotate/eggnog/

wget \
    -nH \
    --user-agent=Mozilla/5.0 \
    --relative \
    --no-parent \
    --reject "index.html*" \
    --cut-dirs=4 \
    -e robots=off \
    -O test-data/funannotate/eggnog/eggnog.db.gz \
    http://eggnogdb.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
    
gunzip test-data/funannotate/eggnog/eggnog.db.gz


wget \
    -nH \
    --user-agent=Mozilla/5.0 \
    --relative \
    --no-parent \
    --reject "index.html*" \
    --cut-dirs=4 \
    -e robots=off \
    -O test-data/funannotate/eggnog/eggnog_proteins.dmnd.gz \
    http://eggnogdb.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz

gunzip test-data/funannotate/eggnog/eggnog_proteins.dmnd.gz
```

## set up interproscan

The easiest is to build a new interproscan container with the DB included,
which takes about 10 hours. See
<https://github.com/TomHarrop/container-interproscan>
