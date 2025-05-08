# Convectively-Coupled Systems in the Tropics as Simulated in Global Storm Resolving Models <img src='https://21centuryweather.org.au/wp-content/uploads/Hackathon-Image-WCRP-Positive-1536x736.jpg' align="right" height="139"/>

Convectively-coupled large-scale systems in the tropics can be categorised into: 1) slow moving moisture modes, 2) fast moving inertio-gravity waves, and 3) mixed systems in between. This is the _**moisture mode-to-gravity wave spectrum**_ framework (see [Adames et al 2019](https://doi.org/10.1175/JAS-D-19-0121.1) and [Adames 2022](https://doi.org/10.1175/JAS-D-21-0215.1)). Moisture mode controls the convection through moistening the tropospheric column, while inertio-gravity wave alters the low-level buoyancy to trigger convection. This project is aimed to find out whether this distinction in how large-scale tropical systems govern the convection is well simulated by the km-scale models.

More details can be found [here](https://github.com/21centuryweather/hk25-teams/blob/main/hk25-AusNode/hk25-AusNode-ConvTrop.md).

**Project leads:**

Martin Singh, Monash University

[Reyhan Respati](mailto:reyhan.respati@monash.edu), Monash University

**Project members:**

name, affiliation/github username

**Collaborators:**

* [hk25-Cyclones](https://github.com/21centuryweather/hk25-teams/blob/main/hk25-AusNode/hk25-AusNode-Cyclones.md)'s works may be useful as we need to exclude tropical cyclones in the analysis.

**Data:**

* Simulation outputs: rlut, pr, ua, va, wa, zg, ta, hus
* Observation datasets: ERA5 (same variables as above), GPM-IMERG precipitation
* Domain: global tropics (30S-30N)
* Time resolution: 3-hourly if possible

## Contributing Guidelines

> Lorem ipsum dolor sit amet 

A python virtual environment has been set up for working on this project. To use the virtual environment, run the following commands everytime you log into `gadi`.

```
module use /g/data/xp65/public/modules
module load conda/analysis3-25.02
source /scratch/gb02/mr4682/tobac_env/bin/activate
```

To close the virtual environment, run this command.

```
deactivate
```

This virtual environment contains all the necessary packages for this project, including `healpy`, `easygems`, and [`tobac`](https://tobac.readthedocs.io/en/latest/), as well as the fundamental libraries like `numpy` and `xarray`.

### Project organisation

All tasks and activities will be managed through GitHub Issues. While most discussions will take place face-to-face, it is important to document the main ideas and decisions on an issue. Issues will be assigned to one or more people and classified using labels. If you want to work on an issue, comment and make sure is assigned to you to avoid overlapping. If you find a problem in the code or a new task, you can open an issue. 

### How to collaborate

* **Main branch:** We want to keep things simple, if you are working on a notebook alone you can push changes to the main branch. Make sure to 1) only add and ccommit that file and nothing else, 2) pull from the remote repo and 3) push.

* **Working on a branch:** if you want to fix or propose a change to someone else code you will need to create a branch and open a pull request. Make sure you explain your suggestion in the pull request message. **This also applies to collaborators outside the project team.**

### Repository structure

This is how the project should look like but make sure to change the name `template-hackathon-project` to something meaningful. 

```bash
template-hackathon-project/
├── LICENCE
├── README.md
├── template_project_hackathon
│   ├── analysis.py
│   ├── __init__.py
│   └── read.py
└── tests
    ├── test_analysis.py
    └── test_read.py
```
* `template_hackathon_project/` this folder will include the code to analysis the data.
* `tests/` this folder contains test code that verifies that your code does what it should.

