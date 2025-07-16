# D2FCSquared

This repository contains MATLAB code and experimental data used to simulate and analyze models described in the accompanying publication. The models include:

- **D2FC²**
- **D2FC optimized**
- **D2FC**
- **IκBβ model**
- **IκBε model**

The experimental data includes:
- The average IKK and nuclear RelA trajectories used for particle swarm optimization in figure 4 and associated validation analysis.  

## Requirements

- **MATLAB**: Version R2024b  
- **SimBiology Toolbox**: Version 24.2 or later  

## Installation  
To install extract to Matlab working folder and navigate to the folder within Matlab. The install time should take less than a minute. 

## Getting Started

To simulate nuclear RelA translocation dynamics based on the **average IKK trajectories**, open the file:

```
RunAverageIKKTrajectory.m
```

in MATLAB and click **Run**. The default model is **D2FC²**. You can change the model by editing the `modelType` variable as explained in the comments at the top of the script. Using the **D2FC²** the expected outcome is a pair of figures displaying the experimental and simulation results of the fitting and validation conditions. Expected runtime is listed below. 

To simulate **single-cell trajectories**, run:

```
RunSingleCellTrajectories.m
```

This also defaults to the **D2FC²** model. To change the model, update the `modelType` variable as described above. Using the **D2FC²** the expected outcome is figures displaying the simulation results of all single cell conditions as in figure 5a of the published article. Expected runtime is listed below. 

All necessary data files for running these simulations are included in this repository.

## Performance

The scripts were tested on a **MacBook Air (Apple M1, 8GB RAM)** and is expected to run simaiarly on other modern laptop and desktop computers:

- `RunAverageIKKTrajectory.m`: ~6.4 seconds  
- `RunSingleCellTrajectories.m`: ~48.9 seconds

## Citation 
Please cite the updated version of: 
"Time-varying stimuli that prolong IKK activation promote nuclear remodeling and mechanistic switching of NF-κB dynamics", Nature Communications, In Press 
