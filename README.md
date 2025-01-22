# TalinStretching
Simulation codes used to investigate work done by a single Talin molecule as it is stretched, its domains unfold, and it unbinds. <br>
(Yamashiro S., Rutkowski D.M., Lynch K.A., Liu Y., Vavylonis D., and Watanabe, N. <i>Force transmission by retrograde actin flow-induced dynamic molecular stretching of Talin.</i> Nature Communications (2023) https://doi.org/10.1038/s41467-023-44018-z). <br>
Talin is represented as beads connected by springs that can unfold by force, bound at one end to a static integrin and the other to a moving actin layer. 

Compilation of pulling simulations can be achieved by using the GNU compiler from within the PullingSimulation directory as: g++ -fopenmp -O3 -o main *.cpp
