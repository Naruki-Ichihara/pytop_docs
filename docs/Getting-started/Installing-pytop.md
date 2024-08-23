---
sidebar_position: 1
---

# Installing pytop

:::warning
pytop do not support python in windows. Please use WSL2.
:::

## Dependencies
This software requires following software.

* [FEniCS 2019](https://fenicsproject.org/download/archive/) (For finite element analysis)
* [NLOpt](https://nlopt.readthedocs.io/en/latest/#download-and-installation) (For optimization)

Abovementioned software cannot be installed using simple ```pip``` command. You need to follow instructions of each softwares.

## On docker contaner
Altanatively, to use our prebuilt container image, you first install [docker CE](https://www.docker.com/products/docker-desktop/), and clone pytop.
```bash
git clone https://github.com/Naruki-Ichihara/pytop.git
```
Then, run ```docker compose```.
```bash
docker compose up
```
Our container will be pulled automatically and running.
:::warning
The container size is over 10GB.
:::
Use *Devcontainer* of [VS code](https://code.visualstudio.com/insiders/) to attach the runnning container. 
Finally, you need to install ```pytop``` in the container, with following command.
 ```bash
pip install .
```

## On codespaces
[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/Naruki-Ichihara/pytop)

Github condespaces is cloud-based virtual environment. Above button navigates and construct virtual vscode with our container on your browser.