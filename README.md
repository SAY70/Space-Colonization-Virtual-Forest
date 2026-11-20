# 🌳 Procedural Tree Generation using MATLAB

## 📚 Overview
This MATLAB script implements a **procedural tree generation** algorithm that creates a forest of trees with varying structural properties. It employs a stochastic approach to generate tree parameters, simulate branching structures, and visualize the results in 3D. The model is designed for use in **radiorealistic vegetation modeling** and **LiDAR-based simulations**.

## 🛠 Key Features
- **Tree Population Control**: Adjust the number of trees, height range, and bud distribution.
- **Branching System**: Categorizes branches into **trunk, primary, secondary, tertiary**, and **non-use branches** based on their radius.
- **Randomized Positioning**: Uses a structured grid with offsets to distribute trees in a defined area.
- **3D Visualization**: Plots tree structures in a 3D environment with labeled receivers.
- **Receiver Perspective Analysis**: Simulates a receiver’s view of the forest and identifies visible branches.

## 🎯 Research Application
This script can be utilized in **remote sensing, forest ecology, and radiative transfer modeling** to analyze the impact of tree structures on wave propagation and environmental factors.

## 🔄 Workflow
1. **Initialize Parameters**: Define tree count, branch magnitude, height limits, and spatial distribution.
2. **Generate Trees**: Use randomized methods to define branching structures based on biological models.
3. **Classify Branches**: Segregate tree components based on diameter thresholds.
4. **Visualize Forest**: Plot trees in 3D with color-coded branches and receiver positions.
5. **Receiver Perspective Simulation**: Identify which tree components fall within the receiver’s cylindrical viewing range.

## ▶️ Getting Started

### Files
| File | Purpose |
|------|---------|
| `spcecol_randdist_tree.m` | Main script: sets the parameters, grows the forest, classifies the branches, plots the 3D scene and runs the receiver analysis |
| `SC_tree_obj.m` | `SC_tree_obj` class used by the script to grow each tree (space colonization) and compute its branch geometry and properties |

### Run
1. Clone or download this repository and keep both `.m` files in the same folder (or add the folder to the MATLAB path).
2. In MATLAB, change to that folder and run:
   ```matlab
   spcecol_randdist_tree
   ```
3. To change the forest, edit the parameters at the top of `spcecol_randdist_tree.m` (see the table below).

> **Note:** the script starts with `clearvars; close all; clc`, so it clears your workspace and closes open figures. It also fixes the random seed with `rng(123774)`, so repeated runs give the same forest until you change the seed.

## 🌟 Notable Parameters
| Parameter | Description |
|-----------|-------------|
| `num_trees` | Total number of trees in the scene |
| `min_h`, `max_h` | Minimum and maximum tree height |
| `min_buds`, `max_buds` | Number of buds per tree |
| `cone_r`, `mid_r`, `bottom_r` | Defines the tapering of tree cones |
| `num_cells_x`, `num_cells_y` | Defines spatial grid for tree distribution |
| `receiver_positions_manual` | Predefined receiver positions in the 3D space |

## 👨‍👩‍👦 Contributing
- Research & Development
- Algorithm Optimization & Testing

## 🔧 Future Enhancements
- **Integration with LiDAR data** for real-world calibration
- **Deep learning models** for tree classification
- **Radiative transfer analysis** for canopy interaction scattering studies

---
🎨 *"Nature is an engineer of complexity. Let's simulate it."* 🌳

