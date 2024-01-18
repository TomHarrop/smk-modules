#!/bin/bash

# kdir -p $HOME/.interproscan-5 && sed 's|^\(data.directory=\).*$|\1/mnt/custom-indices/interproscan/5.59-91.0/data|' $(dirname $(readlink -f $(command -v interproscan.sh)))/interproscan.properties > $HOME/.interproscan-5/interproscan.properties && export _JAVA_OPTIONS=-Duser.home=$HOME &&  interproscan.sh  -dp --input '/mnt/data08/f/f/b/dataset_ffbcff53-5e07-4b99-a5eb-8d18f50

java_home=$1
printf "java_home: %s\n" "${java_home}"
export _JAVA_OPTIONS=-Duser.home=${java_home}

prop_file="$(dirname $(readlink -f $(command -v interproscan.sh)))/interproscan.properties"
printf "prop_file: %s\n" "${prop_file}"



