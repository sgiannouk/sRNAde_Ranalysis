### Maintainer Stavros Giannoukakos 
## To run it: sudo python3 r_installation_for_de.py

## When everything is done, you should be able to access R or Rscript by typing: 
## R_local or Rscript_local

import subprocess, sys, os

r_version = "3.5.1"
main_folder = "/opt/local"

# Creating the directory where R is going to be installed
if not os.path.exists(main_folder): os.makedirs(main_folder)
subprocess.run("sudo chmod 775 -R {0}".format(main_folder), shell=True)

# Installing all necessary dependencies of R
print(">>> Installing all the necessary dependencies of R...")
subprocess.run("sudo apt-get install libpq-dev", shell=True)
subprocess.run("sudo apt-get install libssl-dev", shell=True)
subprocess.run("sudo apt-get install libssh2-1-dev", shell=True)
subprocess.run("sudo apt-get install zlib1g-dev", shell=True)
subprocess.run("sudo apt-get install libcurl4-openssl-dev", shell=True)
subprocess.run("sudo apt-get install phantomjs", shell=True)


# Installing R
print(">>> Installation R version {0}\nThis step might take a while...".format(r_version))
os.chdir(main_folder)
subprocess.run("sudo wget https://cran.rstudio.com/src/base/R-3/R-{0}.tar.gz".format(r_version), shell=True)
subprocess.run("sudo tar xvf R-{0}.tar.gz".format(r_version), shell=True)
os.remove(os.path.join(main_folder, "R-{0}.tar.gz".format(r_version)))
os.chdir(os.path.join(main_folder, "R-{0}".format(r_version)))
subprocess.run("sudo ./configure --prefix={0}/R-{1} --enable-R-shlib --with-blas --with-lapack".format(main_folder, r_version), shell=True)
subprocess.run("sudo make", shell=True)
subprocess.run("sudo make install", shell=True)

# ### R_LIBS_USER=/opt/local/R-3.5.1/library R

# Installing necessary libraries
print(">>> Installation of the necessary R libraries\nThis step might also take a while...")
# print(">> Installation of packrat version 0.5.0")
# subprocess.run("{0}/R-{1}/bin/Rscript -e \'install.packages(\"packrat\", version = \"0.5.0\", lib=\"{0}/R-{1}/library\", repos=\"http://cran.uk.r-project.org\")\'".format(main_folder, r_version), shell=True)
print(">> Installation of ggplot2 version 3.1.0")
subprocess.run("sudo {0}/R-{1}/bin/Rscript -e \'install.packages(\"ggplot2\", version = \"3.1.0\", lib=\"{0}/R-{1}/library\", repos=\"http://cran.uk.r-project.org\")\'".format(main_folder, r_version), shell=True)
print(">> Installation of RColorBrewer version 1.1-2")
subprocess.run("sudo {0}/R-{1}/bin/Rscript -e \'install.packages(\"RColorBrewer\", version = \"1.1-2\", lib=\"{0}/R-{1}/library\", repos=\"http://cran.uk.r-project.org\")\'".format(main_folder, r_version), shell=True)
print(">> Installation of viridis version 0.5.1")
subprocess.run("sudo {0}/R-{1}/bin/Rscript -e \'install.packages(\"viridis\", version = \"0.5.1\", lib=\"{0}/R-{1}/library\", repos=\"http://cran.uk.r-project.org\")\'".format(main_folder, r_version), shell=True)
print(">> Installation of heatmaply version 0.15.2")
subprocess.run("sudo {0}/R-{1}/bin/Rscript -e \'install.packages(\"heatmaply\", version = \"0.15.2\", lib=\"{0}/R-{1}/library\", repos=\"http://cran.uk.r-project.org\")\'".format(main_folder, r_version), shell=True)
print(">> Installation of plotly version 4.8.0")
subprocess.run("sudo {0}/R-{1}/bin/Rscript -e \'install.packages(\"plotly\", version = \"4.8.0\", lib=\"{0}/R-{1}/library\", repos=\"http://cran.uk.r-project.org\")\'".format(main_folder, r_version), shell=True)
print(">> Installation of BiocManager version 1.30.4")
subprocess.run("sudo {0}/R-{1}/bin/Rscript -e \'install.packages(\"BiocManager\", version = \"1.30.4\", lib=\"{0}/R-{1}/library\", repos=\"http://cran.uk.r-project.org\")\'".format(main_folder, r_version), shell=True)
print(">> Installation of vsn version 3.50.0")
subprocess.run("sudo {0}/R-{1}/bin/Rscript -e \'paste(library(\"BiocManager\", lib=\"{0}/R-{1}/library\"), BiocManager::install(\"vsn\", version = \"3.8\", ref=\"3.50.0\", update=FALSE, lib=\"{0}/R-{1}/library\"))\'".format(main_folder, r_version), shell=True)
print(">> Installation of DESeq version 1.34.0")
subprocess.run("sudo {0}/R-{1}/bin/Rscript -e \'paste(library(\"BiocManager\", lib=\"{0}/R-{1}/library\"), BiocManager::install(\"DESeq\", version = \"3.8\", ref=\"1.34.0\", update=FALSE, lib=\"{0}/R-{1}/library\"))\'".format(main_folder, r_version), shell=True)
print(">> Installation of DESeq2 version 1.22.1")
subprocess.run("sudo {0}/R-{1}/bin/Rscript -e \'paste(library(\"BiocManager\", lib=\"{0}/R-{1}/library\"), BiocManager::install(\"DESeq2\", version = \"3.8\", ref=\"1.22.1\", update=FALSE, lib=\"{0}/R-{1}/library\"))\'".format(main_folder, r_version), shell=True)
print(">> Installation of edgeR version 3.24.1")
subprocess.run("sudo {0}/R-{1}/bin/Rscript -e \'paste(library(\"BiocManager\", lib=\"{0}/R-{1}/library\"), BiocManager::install(\"edgeR\", version = \"3.8\", ref=\"3.24.1\", update=FALSE, lib=\"{0}/R-{1}/library\"))\'".format(main_folder, r_version), shell=True)
print(">> Installation of NOISeq version 2.26.0")
subprocess.run("sudo {0}/R-{1}/bin/Rscript -e \'paste(library(\"BiocManager\", lib=\"{0}/R-{1}/library\"), BiocManager::install(\"NOISeq\", version = \"3.8\", ref=\"2.26.0\", update=FALSE, lib=\"{0}/R-{1}/library\"))\'".format(main_folder, r_version), shell=True)

print("\n>>>> The installation has been completed. Please use the following commands in the bash_profiles to create the aliases:")
print("Type: nano ~/.bash_aliases")
print("Copy and past the following two aliases:")
print("alias Rscript_local=\'R_LIBS_USER={0}/R-{1}/library {0}/R-{1}/bin/Rscript\'".format(main_folder, r_version))
print("alias R_local=\'R_LIBS_USER={0}/R-{1}/library {0}/R-{1}/bin/R\'".format(main_folder, r_version))
print("Save and close the file. Then type: source ~/.bash_aliases")

# Exporting the local Rscript and R aliases
# subprocess.run("alias Rscript_local=\'R_LIBS_USER={0}/R-{1}/library {0}/R-{1}/bin/Rscript\'".format(main_folder, r_version), shell=True)
# subprocess.run("alias R_local=\'R_LIBS_USER={0}/R-{1}/library {0}/R-{1}/bin/R\'".format(main_folder, r_version), shell=True)

### Information regarding the local installation of R can be found in the following sources:
# 1. https://support.rstudio.com/hc/en-us/articles/218004217-Building-R-from-source


