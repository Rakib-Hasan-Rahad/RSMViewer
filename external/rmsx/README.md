# RNAMotifScanX Runtime

This folder holds the RNAMotifScanX software that RSMViewer uses when
`data_mode` is `"run_from_scratch"` in `config/rmsx_config.json`. It is **empty
on purpose** in the repository: the first time you run RNAMotifScanX from scratch,
RSMViewer downloads all the required files and third-party software by itself
from

```text
https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx.tar.gz
```

(about 30 MB, verified against `rmsx.tar.gz.sha256`) into this folder and prepares
them for your computer. Later runs start immediately.

Nothing here needs to be installed by hand. Just set `data_mode` to
`run_from_scratch` and run:

```text
rmv_fetch 1S72
rmv_db RNAMotifScanX
```

The default `preannotated` mode does not use this folder at all: it downloads
precomputed annotations per PDB into `external/rmsx_preannotated/`.

`rmv_rmsx_doctor` shows whether everything is ready, and `rmv_setup RNAMotifScanX`
prepares it ahead of time. For requirements on macOS, Windows and Linux, and for
troubleshooting, see [../rmsx_setup.md](../rmsx_setup.md).
