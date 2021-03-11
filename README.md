# James_Foundry_Tools
 
## biomek

#### main.py
	Top level interface for running backdilutions on Biomek, Equal_OD_dilutions
	function takes the output from a plate reader along with a dilution factor,
	target OD, target volume & other variables as inputs and outputs a csv file
	formatted for the biomek's "Transfer from file" step to the "output files"
	folder.
	Prints a warning in case of a well's OD being too low and drops them from
	the experiment  

#### dilution_calculator.py
	Handles calculations related to diluting cultures for main.py

#### biomek_file_writer.py
	general purpose csv writer formatted for the biomek, breaks down transfers
	that exceed maximum allowed volume into multiple smaller transfers
	automatically, can save to .csv or return csv formatted file to console

## flow_cytometry

#### BD_accuri_parser.py
	built off the open source project "FlowIO", reads .fcs files outputted by
	the flow cytometer and combines with a strain map to provide easily
	graphable and retrievable results
	can make multistrain comparing any combination of two measurements, across
	an arbitrary number of strains. Can take custom dimensions for side by side
	comparison or generate dimensions automatically
	Can also retrieve results for a single strain for custom figures.

## microtiter

#### plate_reader_tools.py
	contains an array of basic tools used throughout the package including:
		Reading one-shot plate reader outputs
		Reading growth curve data from sunrises 
		Conversion between matrix format and list format
		Conversion between well coordinates (A1-H12) and numbers (1-96)
		Reading strain & treatment maps for other modules

#### curve_maker.py
	Contains higher level tools for analyzing growth curves.
	curve_maker accepts a strain map path and a sunrise path and outputs a
	long format dataframe with columns as Time, well, name & OD. Ideal for
	combining output of multiple experiments into a single dataframe
	curve_viewer recieves the output of curve_maker (works better with
	raw_time set to False) and produces a line plot using the seaborn.relplot
	function. 
		accepts an optional  context variable with “notepad”, “paper”, “talk”,
		and “poster” as options,
		accepts optional: legend_name, and a list of names from the strain map
		to plot
	mu_max recieves the output from curve_maker and produces the mu max of each
	well, retaining identifier information 
	mu_max_plot provides a high level function for plotting mu max results,
	designed specifically for multiple strains & treatments with a control
	strain 

