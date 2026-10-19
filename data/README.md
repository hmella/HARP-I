## Example datasets

`data_2D.mat` is **already in this folder**, committed with the repository.
The 2D demos and the Python test suite work from a fresh clone without
downloading anything.

`data_3D.mat` is 1.3 GB and is too large to keep in git, so it is hosted
externally:

🔗 **[Download the datasets](https://www.dropbox.com/scl/fo/fs73nh11mpl2krbfytmu6/AHQR7jUvlGHwE1PY2csGBtE?rlkey=e18m8dsmdtfv9anrp90dv2ttk&st=pn7h1jzd&dl=0)**

Copy `data_3D.mat` into this folder. Without it the 3D demos cannot run and a
few Python tests skip, naming the missing file; nothing else is affected.

### Tagged DICOM examples

Some Python tests read clinical DICOM exports from a `Tagging dataset/`
folder here. These are patient images and are not distributed. Those 37 tests
skip when the folder is absent, which is the normal case, and `harpi dicom
<path>` reads any export of your own.
