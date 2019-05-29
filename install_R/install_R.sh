#!/usr/bin/env bash

if [[ $EUID -ne 0 ]]; then
   echo "This script must be run as root, use sudo "$0" instead" 1>&2
   exit 1
fi

r_url="https://cran.rstudio.com/src/base/R-3/R-3.5.3.tar.gz"
current_path=$(readlink -f $0)
current_dir=$(dirname $current_path)
main_dir="/opt/local"


if [ ! -d "$main_dir" ]; then
    mkdir -p "$main_dir"
fi

apt update && apt install -y wget curl build-essential gfortran libcurl4-openssl-dev libxml2-dev libbz2-dev liblzma-dev libpcre3 libpcre3-dev libpq-dev libssl-dev libssh2-1-dev zlib1g-dev libcurl4-openssl-dev phantomjs
cd "$main_dir"
wget -c "$r_url" && tar xfv R-3.5.3.tar.gz && rm R-3.5.3.tar.gz
cd R-3.5.3
./configure --enable-R-shlib --with-blas --with-lapack --with-readline=no --with-x=no && make && make install
cd -
"$main_dir/R-3.5.3/bin/Rscript" "$current_dir/install_dependencies.R"
