# Weibull Waiting Time Analysis

This repository contains a Python pipeline for analyzing waiting time distributions of burst events (e.g., Fast Radio Bursts, Giant Pulses) using a Weibull model. It utilizes [Cobaya](https://cobaya.readthedocs.io/) for Bayesian inference to constrain the Weibull shape parameter ($k$) and rate parameter ($r$).

## Features

- **Weibull Modeling**: Fits waiting times between bursts to a Weibull distribution.
- **Bayesian Inference**: Uses MCMC integration creates posterior chains for model parameters ($k$ and $r$). The rate parameter $r$ is measured in $\mathrm{hr}^{-1}$.
- **Flexible Input**: Accepts custom observation windows and burst timing data.
- **Automated Plotting**: Generates corner plots of the posterior distributions using [GetDist](https://getdist.readthedocs.io/).

## Installation

Ensure you have Python 3 installed. Install the required dependencies:

```bash
pip install cobaya scipy numpy getdist matplotlib
```

## Usage

### 1. Prepare Data
The analysis requires a single **JSON** file containing a list of observation epochs. Each element in the list must specify the `start` and `end` times (MJD) of the epoch and a list of `bursts` times (MJD).

**Format**:
```json
[
  {
    "id": 1,
    "start": 59593.709,
    "end": 59593.729,
    "bursts": [59593.71119, 59593.72153]
  },
  {
    "id": 2,
    "start": 59594.100,
    "end": 59594.150,
    "bursts": []
  }
]
```

### 2. Configure
Edit `config.yaml` to point to your data file:

```yaml
likelihood:
  waiting_time_likelihood.WeibullLikelihood:
    data_file: "path/to/your/data.json"

params:
  k:
    prior:
      min: 0.01
      max: 10
...
```

### 3. Run
Execute the runner script:

```bash
python run_cobaya.py
```

### 4. Output
Results are saved in the `chains/` directory:
- MCMC chains (`.txt` files).
- Corner plot PDF (e.g., `chains/weibull_corner.pdf`) showing the marginalized constraints on $k$ and $\log_{10}r$.

## Project Structure

- `run_cobaya.py`: Entry point script. Runs the sampler and plotting.
- `config.yaml`: Configuration file for Cobaya (inputs, priors, sampler settings).
- `waiting_time_likelihood.py`: Custom Cobaya likelihood module implementing the Weibull model logic.

## Methodological Caveats

This analysis is based on the methods described in [Oppermann, Yu, & Pen (2018)](https://academic.oup.com/mnras/article/475/4/5109/4791594). Users should be aware of the following assumptions and limitations:

1.  **Independence of Epochs**: The analysis treats each observation epoch as statistically independent. For a Weibull process with $k \neq 1$ (non-Poissonian), the process is not memoryless. Treating epochs as independent neglects the "memory" of the process (time since the last burst) carried over the unobserved gaps between sessions. The likelihood calculation for the first burst in an observation assumes the process "resets" or uses a standard censoring approximation at the start of the window.
2.  **Stationarity**: The model assumes the underlying physical parameters ($k$ and $r$) are constant over the entire timespan of the dataset. It models the "bursty" nature (active vs. quiescent phases) via the clustering parameter $k$ rather than explicit on/off states.
3.  **Sensitivity Dependence**: The derived parameters ($r$ and $k$) are effective values for the specific sensitivity threshold of the observations. If sensitivity varies significantly between epochs or if many weak bursts are missed, the observed waiting time distribution is a convolution of the intrinsic distribution and the detection efficiency.

