# Instructions

## Conda packages required

```bash
conda install anaconda-client conda-build
```

## Building and pushing to https://anaconda.org/uibcdf

```bash
conda build . --no-anaconda-upload --python 3.7
PACKAGE_OUTPUT=`conda build . --output`
anaconda login
anaconda upload --user uibcdf $PACKAGE_OUTPUT
conda build purge
anaconda logout
```
## Install

```
conda install -c uibcdf openmolecularsystems
```

## Additional Info
https://docs.anaconda.com/anaconda-cloud/user-guide/tasks/work-with-packages
