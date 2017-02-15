# Usr directory 
Packages are installed in the usr directory. 

## Packages 

### MaBoSS - boolean network simulator
MaBoSS can be run on Linux or Windows OS, but I will explain the installation process for Linux OS. You can download, build, and install MaBoSS by running the following commands:

```bash 
wget https://maboss.curie.fr/pub/MaBoSS-env-2.0.tgz
cd MaBoSS-env-2.0/
./check-requirements.sh 
cd engine/src
make 
# or 
make install
# Then, you can find the bin file MaBoSS.  
```

After installation, run the following commands to enable MaBoSS to run.

```bash 
cd MaBoSS-env-2.0/
source MaBoSS.env 
```

If you do not want to run the above command for each time, add the following code lines to .bashrc:

```bash 
TOPDIR=/home/pbs/git/sbie_aging/usr/MaBoSS-env-2.0
tooldir=tools
PATH=$TOPDIR/${tooldir}:$TOPDIR/engine/pub:$PATH; export PATH
PERL5LIB=$TOPDIR/${tooldir}:${PERL5LIB}; export PERL5LIB
PYTHONPATH=$TOPDIR/${tooldir}:${PYTHONPATH}; export PYTHONPATH
```


