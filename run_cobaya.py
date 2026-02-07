from cobaya.run import run
from cobaya.yaml import yaml_load_file
import os
from getdist import plots, loadMCSamples

if __name__ == "__main__":
    if not os.path.exists("chains"):
        os.makedirs("chains")

    info = yaml_load_file("config.yaml")
    updated_info, sampler = run(info)

    # Plotting
    output_prefix = updated_info.get("output")
    if output_prefix:
        print(f"Generating plot for chains: {output_prefix}")
        try:
            samples = loadMCSamples(output_prefix, settings={'ignore_rows': 0.3})
            g = plots.get_subplot_plotter()
            g.triangle_plot(samples, ["k", "logr"], filled=True, title_limit=1)
            
            plot_filename = output_prefix + "_corner.pdf"
            g.export(plot_filename)
            print(f"Corner plot saved to {plot_filename}")
        except Exception as e:
            print(f"Error generating plot: {e}")
