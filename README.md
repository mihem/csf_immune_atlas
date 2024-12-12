# CSF Immune Atlas Code
This repository contains the code used in the manuscript
[Atlas of cerebrospinal fluid immune cells across neurological diseases](http://doi.org/10.1002/ana.27157)
published in *Annals of Neurology* 2024.

# Structure

| Directory            | Description |
| ---------            | -----------  |
| `scripts/analysis`   | Code used for the analysis in the manuscript |
| `models`             | Final XGB models |

# Package for helper functions
Helper functions of this analysis are bundled in the package [CSFAtlasTools](https://github.com/mihem/CSFAtlasTools)

# Reproducibility
To ensure reproducibility, we used the *renv* package and Docker. To restore the packages from the *renv.lock* file, use:

```R
renv::restore()
```
Alternatively, you can also use the Docker image, which contain all necessary R packages and system dependencies. The Docker image is available on DockerHub:

```bash
docker pull mihem/csf_immune_atlas:2.0
```


# Contact
If you have any questions, please contact me via [mheming.de](https://osmzhlab.uni-muenster.de/mheming/#contact).