# Pixels and their neighbours: Finite volume

_Rowan Cockett, Lindsey Heagy and Doug Oldenburg_ - [The Leading Edge](http://library.seg.org/doi/abs/10.1190/tle35080703.1)

[![Made with MyST](https://img.shields.io/badge/made%20with-myst-orange)](https://curve.space/examples/pixels)
[![Binder](http://mybinder.org/badge.svg)](http://mybinder.org:/repo/simpeg/tle-finitevolume)
[![DOI](https://img.shields.io/badge/doi-10.1190%2Ftle35080703.1-blue)](https://doi.org/10.1190/tle35080703.1)

These notebooks use [Python 3.8](https://docs.python.org/) and the open source packages [SimPEG](https://simpeg.xyz) and [discretize](https://discretize.simpeg.xyz).
[SimPEG](https://simpeg.xyz) can be installed using the Python package manager `PyPi` or `conda` and running:

```
conda install SimPEG --channel conda-forge
```

**Overview**

In [index.ipynb](/notebooks/index.ipynb), we first give an overview and introduce the problem, in [mesh.ipynb](/notebooks/mesh.ipynb) we show how to create and work with a regular mesh, in [divergence.ipynb](/notebooks/divergence.ipynb) we build the discrete divergence operator, and in [weak-formulation.ipynb](/notebooks/weak-formulation.ipynb), we discretize and solve the DC equations using weak formulation. Finally we bring it all together with an interactive app in [all-together-now.ipynb](/notebooks/all-together-now.ipynb)

**Repository Contents**

- [Paper](paper.md)
- Notebooks:
  - [index](/notebooks/index.ipynb)
  - [mesh](/notebooks/mesh.ipynb)
  - [divergence](/notebooks/divergence.ipynb)
  - [weak-formulation](/notebooks/weak-formulation.ipynb)
  - [all together now](/notebooks/all-together-now.ipynb)

**Cite Us 🙏**

> Rowan Cockett, Lindsey J. Heagy, and Douglas W. Oldenburg (2016).
> "Pixels and their neighbors: Finite volume." The Leading Edge, 35(8), 703–706.
> https://doi.org/10.1190/tle35080703.1

**Made With MyST**

This paper and notebooks are made with [MyST Markdown](https://myst.tools). Which creates scientific documents, and can publish to the web (see [on Curvenote](https://next.curve.space/examples/pixels)).

[![](/images/screenshot.png)](https://next.curve.space/examples/pixels)
