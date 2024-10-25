# CSF immune atlas code
This repository contains the code used in the manuscript "Atlas of cerebrospinal fluid immune cells across neurological diseases".
Published in ... DOI

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

To build the Docker image, use:

```bash
docker pull mihem/csf_immune_atlas:1.0
```


# Contact
If you have any questions, please contact me via [mheming.de](https://osmzhlab.uni-muenster.de/mheming/#contact).