# Installation

## Latest stable version

The latest stable version of UIBCDF_test_systems can be installed from the UIBCDF Anaconda channel:

```bash
conda -c uibcdf uibcdf_test_systems
```

Once the library was installed, it can be removed from your conda environment with:

```bash
conda remove uibcdf_test_systems
```

## Developing version

The public repository of the whole library is available in the UIBCDF GitHub page. It can be cloned
in your computer to install it locally (within or without a conda environment, in this case is
indifferent).

```bash
git clone https://github.com/uibcdf/UIBCDF_test_systems.git
cd uibcdf_test_systems
python setup.py develop
```

To uninstall this version of the library:

```bash
pip uninstall uibcdf_test_systems
```

