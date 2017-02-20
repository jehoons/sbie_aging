# Aging research
This project is to study aging process and how to overcome or reverse it. See also [Google Drive][Google Drive].

## Learning git
To participate in this project, first of all, you need to learn about git. First, install [git-scm][git-scm]. After installation, you can right-click and select `git bash here`. If selected, `bash shell` will be executed at that location. For more information, see also [git-cheat-sheet][git-cheat-sheet]. To make local clone (or repository) of project to your computer, run the following command:

```
# in case that you are using Windows OS
$git clone https://github.com/jehoons/sbie_bone.git

or

# in case that you are using Linux OS
$git clone git@github.com:jehoons/sbie_bone.git
```

## Packages
### MaBoSS - Markovian Boolean Stochastic Simulator

[MaBoSS][MaBoSS-site] can be run on Linux or Windows OS, but I will explain the installation process for Linux OS. You can download, build, and install MaBoSS by running the following commands:

```bash
cd usr
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

If you do not want to run the above command for each time, add the following code lines to `.bashrc`:

```bash
TOPDIR=/home/pbs/git/sbie_aging/usr/MaBoSS-env-2.0
tooldir=tools
PATH=$TOPDIR/${tooldir}:$TOPDIR/engine/pub:$PATH; export PATH
PERL5LIB=$TOPDIR/${tooldir}:${PERL5LIB}; export PERL5LIB
PYTHONPATH=$TOPDIR/${tooldir}:${PYTHONPATH}; export PYTHONPATH
```

See also [format][maboss-format], and [paper][maboss-paper].

[Google Drive]: https://drive.google.com/open?id=0B2Fh-6_aEya5MU9nTldLN2FIVW8
[git-scm]: https://git-scm.com/download/win
[git-cheat-sheet]: https://www.git-tower.com/blog/git-cheat-sheet
[maboss-format]: https://maboss.curie.fr/pub/DescriptionOutputFile.pdf
[maboss-paper]: http://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-6-116
[MaBoSS-site]: https://maboss.curie.fr
