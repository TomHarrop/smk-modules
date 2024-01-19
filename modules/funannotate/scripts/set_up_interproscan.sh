#!/bin/bash

# this is how galaxy runs it:
# mkdir -p $HOME/.interproscan-5 && sed 's|^\(data.directory=\).*$|\1/mnt/custom-indices/interproscan/5.59-91.0/data|' $(dirname $(readlink -f $(command -v interproscan.sh)))/interproscan.properties > $HOME/.interproscan-5/interproscan.properties && export _JAVA_OPTIONS=-Duser.home=$HOME &&  interproscan.sh  -dp --input '/mnt/data08/f/f/b/dataset_ffbcff53-5e07-4b99-a5eb-8d18f5018610.dat' --seqtype p -f TSV,XML  --applications TIGRFAM,FunFam,SFLD,SUPERFAMILY,PANTHER,Gene3D,Hamap,PrositeProfiles,Coils,SMART,CDD,PRINTS,PIRSR,PrositePatterns,AntiFam,Pfam,MobiDBLite,PIRSF --tempdir ${TEMP:-$_GALAXY_JOB_TMP_DIR}  --pathways --goterms   --cpu ${GALAXY_SLOTS:-4}  --output-file-base 'output'

java_home=$1
printf "java_home: %s\n" "${java_home}"
export _JAVA_OPTIONS=-Duser.home=${java_home}

prop_file="$(dirname $(readlink -f $(command -v interproscan.sh)))/interproscan.properties"
printf "prop_file: %s\n" "${prop_file}"



