from mono_analysis import plot_raw_traces
import gcms_analysis

if __name__ == '__main__':

	plot_raw_traces.make_species_di_mono_carb_plot()
	gcms_analysis.compile_gcms_stats_to_pandas.make_solution_size_summary_figure()

	pd = gcms_analysis.compile_gcms_stats_to_pandas.compile_all_data()
	gcms_analysis.compile_gcms_stats_to_pandas.make_correlation_matrix_figure(pd)
	gcms_analysis.compile_gcms_stats_to_pandas.update_numerical_chc_values(pd)
	#gcms_analysis.compile_gcms_stats_to_pandas.plot_correlation_with_size_removed(pd, 'carb', ['hairyness', 'RetTimeAvg', 'pulvilli'])

	plot_raw_traces.make_frozen_hexane_plot()
	