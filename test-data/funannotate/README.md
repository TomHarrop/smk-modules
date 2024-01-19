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
    docker://ghcr.io/tomharrop/container-funannotate:1.8.15_cv2 \
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
    docker://ghcr.io/tomharrop/container-funannotate:1.8.15_cv2 \
        funannotate setup \
            -i all \
            -d test-data/funannotate/db
```

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

## get the eggnog DB from the container

```bash
tmpdir="$( mktemp -d )"
IPR_DIR="$( readlink -f test-data/funannotate/interproscan )"

export JAVA_OPTIONS=-Duser.home=$IPR_DIR

mkdir -p "${IPR_DIR}"
mkdir -p "${IPR_DIR}/.interproscan-5"

apptainer exec \
    -B ${tmpdir} \
    -H ${tmpdir} \
    --pwd ${tmpdir} \
    --containall \
    --cleanenv \
    --writable-tmpfs \
    docker://quay.io/biocontainers/interproscan:5.59_91.0--hec16e2b_1 \
        cp -r /usr/local/share/InterProScan/data \
        /usr/local/share/InterProScan/interproscan.properties \
        ./ 

mv "${tmpdir}/data" "${IPR_DIR}/"

sed \
"s|^\(data.directory=\).*$|\1${IPR_DIR}/data|" \
"${tmpdir}/interproscan.properties" \
> "${IPR_DIR}/.interproscan-5/interproscan.properties"

# Run hmmpress on all hmm files in the interproscan directory.
# This is what setup.py is supposed to do.
apptainer exec \
    -B ${IPR_DIR} \
    -H ${tmpdir} \
    --pwd ${tmpdir} \
    --containall \
    --cleanenv \
    --writable-tmpfs \
    --env IPR_DIR=$IPR_DIR \
    docker://quay.io/biocontainers/interproscan:5.59_91.0--hec16e2b_1 \
    find ${IPR_DIR}/data -type f -name "*.hmm" \
    -exec /usr/local/bin/hmmpress {}  \; 

# Test run interproscan
    --containall \
    --cleanenv \
    --writable-tmpfs \


apptainer exec \
    -B ${IPR_DIR} \
    -B ${PWD} \
    -H ${tmpdir} \
    --pwd ${PWD} \
    --env IPR_DIR=$IPR_DIR \
    docker://quay.io/biocontainers/interproscan:5.59_91.0--hec16e2b_1 \
    bash -c '\
        export _JAVA_OPTIONS=-Duser.home=${IPR_DIR} &&
        interproscan.sh \
        -dp \
        -i test-output/funannotate/funannotate/predict_results/testspecies.proteins.fa \
        --output-dir test_ipr \
        --cpu 11 \
        ' > log.out 2> log.err

```
