# physarum_maze
Code to reproduce the numerical results for "Threshold sensing yields optimal path formation in Physarum polycephalum"

## Code:
- The code contained in folders `lib`, `src` and `test` is adapted from a previously optimised ODE solver for similar problems, by David Palma (as acknowledged in the source code). It is hosted for reproducibility purposes of the specific paper results, with tuned parameters and conditions.
- The file _physarum_maze_sim_stability_markers.m_ also contains code to simulate the dynamics in a more complex maze, constructed via recursion. Although not used in the paper, it provides a nice playground to test configurations and build intuition about the mould behaviour.

## Credits:
If you reuse either code with these specific configurations, please cite Daniele Proverbio, Giulia Giordano, "Threshold sensing yields optimal path formation in Physarum polycephalum", arXiv preprint arXiv:2507.12347

If you wish to reuse the ODE solver as is, we recommend citing the original publication: 
Franco Blanchini, Daniele Casagrande, Filippo Fabiani, Giulia Giordano, David Palma, and Raffaele Pesenti. A threshold mechanism ensures minimum-path flow in lightning discharge. Scientific Reports, 11, 12 2021. doi: 10.1038/s41598-020-79463-z
