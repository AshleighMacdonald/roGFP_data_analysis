# roGFP_data_analysis
This is an example R script which was used to analyse data collected from bacteria which had been transformed with reduction-oxidation sensitive GFP (roGFP).
roGFP was used to monitor the intracellular redox state of bacterial gene deletion mutants to different concentrations of a treatment (which is hydrogen peroxide in this example script).
This was completed to assess the roles of genes in response to the treatment.

See my thesis for further details 'Ozone-Mediated Control of Food Spoilage and Food-Borne Pathogens.'

The data which was collected from roGFP experiments was analysed as described by van der Heijden and Finlay, 2015. Briefly, background fluorescence signal was corrected for at each treatment condition by subtracting the fluorescence values measured from the pDSK519-expressing strain, from the fluorescence values measured from the pDSK-roGFP-expressing strain, at each time point. Then, the corrected fluorescence values which were measured at 405 nm excitation were divided by the corrected fluorescence values which were measured at 480 nm excitation at the same time point, to provide a ‘405/480 nm ratio’. Then in every individual assay, the data was normalised using 405/480 nm ratio measurements from the same bacterial strain in a maximum oxidised and maximum reduced state. A normalised 405/480 ratio of 1 was made to be equal to the maximum oxidised signal, and 0.1 was made to be equal to the maximum reduced signal.
