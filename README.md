# CellularAutomata-SnowCrystals
This Repository contains a Matlab code to efficiently reproduce the results of the paper [A local cellular model for snow crystal growth](https://www.sciencedirect.com/science/article/pii/S0960077904003741) by Clifford A.Reiter.
All the codes in this repository are Matlab codes.

## Code Versions:
### v_1:
This code generates the crystal structure but uses hexagons to make each cell which is computationally intensive and computers often run out of RAM if large grids are used. However this is the most accurate version.

### v_2:
This code is much faster as it uses pixels to simulate hexagons and is therefore appreciatably fatser than v_1. Eventhough it is faster it doesnt lose much in overall quality of the results.
