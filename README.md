# 2025-aqgeo-pde-magnetics-tutorial

## How to run the notebooks

### Install a Python distribution

In order to run the notebooks you need to install a Python distribution like
[Miniforge][miniforge], [Anaconda][anaconda] or [Miniconda][miniconda].

Follow the instructions provided by any of them to install it in your system.

> [!NOTE]
> In the following steps we'll instruct you to run some commands. They are
> meant to be run in a terminal (for Mac and Linux users), or in an Anaconda
> Prompt or Miniforge Prompt if you are running Windows.

### Download the repository

Then, download this repository as a [zip file][repo-zip], or clone it with `git`:

```bash
git clone https://github.com/ubcgif/2025-aqgeo-pde-magnetics-tutorial
```

### Navigate to the folder

```bash
cd 2025-aqgeo-pde-magnetic-tutorial
```

### Create a new conda environment

In the downloaded folder you'll find an `environment.yml` file that you can use
to create a conda environment by running:

```bash
conda env create -f environment.yml
```

Once the new environment is created and all the packages are installed, we need
to activate the environment with:

```bash
conda activate aqgeo_pdemag
```

## Run Jupyter Notebook

After we activated the environment, let's start a Jupyter notebook, where we can open,
edit and execute the notebooks:

```bash
jupyter notebook
```

I will automatically open a new tab in your browser with the Jupyter notebook, where
you can navigate into the `notebooks` folder and open any of the notebooks
you'll find there.

> [!IMPORTANT]
> It's important not to close the terminal (or Prompt) window while you are
> using the Jupyter notebook. Closing the terminal will quit the Jupyter notebook.

[anaconda]: https://anaconda.org
[miniconda]: https://docs.anaconda.com/miniconda/miniconda-install
[miniforge]: https://github.com/conda-forge/miniforge
[repo-zip]: https://github.com/ubcgif/2025-aqgeo-pde-magnetics-tutorial/archive/refs/heads/main.zip