The FPoLi-GPCR (Fluorescence Polarization Ligand-GPCR binding) folder contains equations and ready to use example files and software to
analyze fluorescence polarization/fluorescence anisotropy based ligand binding experiments.

When using these models or their derivatives please cite:

Laasfeld, T., Tahk, M. J., Allikalt, A., Torp, J., Grätz, L., Kopanchuk, S., & Rinken, A. (2024). Exploring Muscarinic Acetylcholine Receptor Binding Kinetics with Fluorescence Anisotropy. In Muscarinic Receptor: From Structure to Animal Models (pp. 113-151). New York, NY: Springer US.

Veiksina, S., Kopanchuk, S., & Rinken, A. (2014). Budded baculoviruses as a tool for a homogeneous fluorescence anisotropy-based assay of ligand binding to G protein-coupled receptors: The case of melanocortin 4 receptors. Biochimica et Biophysica Acta (BBA)-Biomembranes, 1838(1), 372-381.

NB! Care must be taken when using pzfx and xml files as different versions of GraphPad Prism converge to different, sometimes ambiguous solutions. Use .pzf file for reference of original fits. If your Prism version fails to converge at particular conditions it may have to do with the way rounding of very small numbers is handled. When using the equations at these conditions/parameters care should be taken.

Contents:
SaturationBinding
	-SaturationBinding.pzf                   - A GraphPad Prism v5 version file showing data organization and analysis parameters for the saturation binding experiments.
	-SaturationBinding.pzfx                  - A GraphPad Prism v5 version file showing data organization and analysis parameters for the saturation binding experiments.
	-SaturationBinding.xml                   - A GraphPad Prism xml file showing data organization and analysis parameters for the saturation binding experiments.
	
	-SaturationBindingFP.txt                 - A simple text file with the saturation binding model equations, useful for implementing in software other than GraphPad Prism
	-SaturationBindingFP.m                   - A MATLAB version of the saturation binding model equations.

CompetitionBinding
	-CompetitionBinding.pzf                  - A GraphPad Prism v5 version file showing data organization and analysis parameters for the competition binding experiments.
        -CompetitionBinding.pzfx                 - A GraphPad Prism v5 version file showing data organization and analysis parameters for the competition binding experiments.
	-CompetitionBinding.xml                  - A GraphPad Prism xml file showing data organization and analysis parameters for the competition binding experiments.
	
	-Wang-Kopanchuk_equation.txt             - A simple txt file with the Wang-Kopanchuk competition binding equation, useful for implementing in software other than GraphPad Prism
	-Wang-Kopanchuk_equation.m               - A MATLAB version of the Wang-Kopanchuk competition binding equation. Allows choosing between the regular, less numerically 
						   exact but faster version, a version using symbolic math or automatically adjusting version
	

IQMTools4Aparecium
	-

Competition Binding Simulator
	
	- 
