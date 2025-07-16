# D2FCSquared

This repository contains MATLAB code used to simulate and analyze models described in the accompanying publication. The models include:

- **D2FC²**
- **D2FC optimized**
- **D2FC**
- **IκBβ model**
- **IκBε model**

## Requirements

- **MATLAB**: Version R2024b  
- **SimBiology Toolbox**: Version 24.2 or later  

## Getting Started

To simulate nuclear RelA translocation dynamics based on the **average IKK trajectories**, open the file:

```
RunAverageIKKTrajectory.m
```

in MATLAB and click **Run**. The default model is **D2FC²**. You can change the model by editing the `modelType` variable as explained in the comments at the top of the script.

To simulate **single-cell trajectories**, run:

```
RunSingleCellTrajectories.m
```

This also defaults to the **D2FC²** model. To change the model, update the `modelType` variable as described above.

All necessary data files for running these simulations are included in this repository.

## Performance

The scripts were tested on a **MacBook Air (Apple M1, 8GB RAM)**:

- `RunAverageIKKTrajectory.m`: ~6.4 seconds  
- `RunSingleCellTrajectories.m`: ~48.9 seconds
