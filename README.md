# Cellcano

---

[![PyPI version](https://badge.fury.io/py/Cellcano.svg)](https://badge.fury.io/py/Cellcano) [![Python 3.8](https://img.shields.io/badge/python-3.8-blue.svg)](https://www.python.org/downloads/release/python-360/) [![DOI](https://zenodo.org/badge/449052687.svg)](https://zenodo.org/badge/latestdoi/449052687)

Cellcano is an open source software for supervised cell type identification (celltyping) in scATAC-seq data, published in [Nature Communications](https://doi.org/10.1038/s41467-023-37439-3). The motivation to develop Cellcano are:
1. Supervised methods are more accurate, robust and efficient than unsupervised clustering methods in scATAC-seq data
2. With more high-quality scATAC-seq datasets being generated, methods using scATAC-seq as references can have better prediction performances and are in high demand

More details and tutorial: https://marvinquiet.github.io/Cellcano/.

**Table of Contents**
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [License](#license)
- [Citing Our Work](#citing-our-work)


## System Requirements

### Hardware requirements

Cellcano package requires only a standard computer with enough RAM to support the in-memory operations. Cellcano can use GPU if the computer has the GPU resource but it is not required.

### Software requirements

#### OS requirements

Cellcano supports `macOS`, `Linux` and `Windows`. It has been tested on all three systems.

#### Dependencies

Cellcano requires the following:

- [python](https://www.python.org/) (3.8 recommended)
- [R](https://www.r-project.org/)
- [tensorflow](https://www.tensorflow.org/) (2.7.1)
- [anndata](https://anndata.readthedocs.io/en/latest/) (0.7.4)
- [scanpy](https://scanpy.readthedocs.io/en/stable/) (1.8.2)
- [numpy](https://numpy.org/) (1.19.2)
- [h5py](https://www.h5py.org/) (2.10.0)
- [keras](https://keras.io/) (version compatible with tensor flow)
- [rpy2](https://rpy2.github.io/) (version compatible with both Python and R)
- `cuda toolkit` and `nvidia cudnn` if using GPU, more information can be found [here](https://towardsdatascience.com/setting-up-tensorflow-gpu-with-cuda-and-anaconda-onwindows-2ee9c39b5c44)

If the input is scATAC-seq raw data (i.e. fragment file or bam file), [ArchR](https://www.archrproject.com/) package has to be installed. 



## Installation

The most convinient way is to install with `pip`.

```shell
pip install Cellcano
```

To upgrade to a newer release use the `--upgrade` flag.

```shell
pip install --upgrade Cellcano
```

We have a detailed tutorial on installation in our [documentation](https://marvinquiet.github.io/Cellcano/).



## License

This project is covered under the **MIT license**.



## Citing Our Work

For usage of the package and associated manuscript, please cite: 
```BibTex
@article{ma23cellcano,
  title   = {Cellcano: supervised cell type identification for single cell ATAC-seq data},
  author  = {Ma, Wenjing and Lu, Jiaying and Wu, Hao},
  journal = {Nature Communications},
  year    = {2023},
  month   = {Apr.},
  day     = {03},
  volume={14},
  number={1},
  pages={1864},
  issn={2041-1723},
  doi={10.1038/s41467-023-37439-3},
  url={https://doi.org/10.1038/s41467-023-37439-3}
}
```

