# Starvation_paper_AT
To share with Jamie, Nathan and Tyler

Similarly named files are modified version of each other.

For AFD imaging, analysis were done..
Ca_response_calc_cell_body -> baseline_correct_universal_AT (on GCaMP_Analysis....mat in Alayisis-files folder created by previous script)-> (pile up blcrrt_Expt...mat files) -> plot_multiple_calcium_heatmap_AT or find_T_star_linear_AT  (on the folder of piled up blcrrt files to make heatmaps or calculate T*) 

For AIY imaging (first and last codes are the only difference from AFD),
Ca_response_calc_neurite (on tiff)-> baseline_correct_universal_AT  -> plot_multiple_calcium_heatmap_AT or Ca_events_GCaMP_base10 (to calculate total duration of response) or imaging_proporion_heatmap (to generate the histogram&get the data for it) (all on the folder of piled up blcrrt files)

For AWC&AIB imaging (Similar to AIY, but cut of the data of temperature range from 20-27ºC)
Ca_response_calc_cell_body -> Ca_response_partial_analysis -> baseline_correct_universal_AT(on partial_GCaMP_Analysis....mat in Alayisis-files folder created by previous script) -> (pile up blcrrt_partial_Expt...mat files) -> plot_multiple_calcium_heatmap_AT_partial or Ca_events_GCaMP_base10_partial or imaging_proporion_heatmap

For AIA&AIZ imaging (Similar to AWC&AIB, but the first script is different)
Ca_response_calc_neurite -> Ca_response_partial_analysis -> baseline_correct_universal_AT(on partial_GCaMP_Analysis....mat in Alayisis-files folder created by previous script) -> (pile up blcrrt_partial_Expt...mat files) -> plot_multiple_calcium_heatmap_AT_partial or Ca_events_GCaMP_base10_partial or imaging_proporion_heatmap
