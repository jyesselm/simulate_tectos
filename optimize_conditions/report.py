import os
import shutil

import plot
import datasets

class Report(object):
    def __init__(self):
        self.fig_count = 0

    def run_plot(self):
        pass

    def get_html(self):
        pass

class TectoBasicReport(Report):
    def __init__(self, dir):
        self.dir = dir

    def get_html(self, df):
        pnames = """TectoOverallCompare,
                    TectoBasepairStepCorrelation,""".split(",")
        pnames = [pname.rstrip().lstrip() for pname in pnames]
        plotters = [plot.PlotFactory.get_plot(pname) for pname in pnames]
        dts = [df]
        dts.append(datasets.basepair_step_dependence(df))
        #dts.append(datasets.basepair_step_position_dependence(df))

        f = open(self.dir+"/report.html", "w")
        f.write("<html>\n")
        f.write("<body>\n")

        for i in range(len(plotters)):
            section_dir = self.dir+"/section."+str(i)
            if os.path.exists(section_dir):
                shutil.rmtree(section_dir)
            os.mkdir(section_dir)
            plotters[i].savefig(dts[i], dir=section_dir)

            f.write("ploting class: " + plotters[i].name + "<br>")
            for j in range(plotters[i].fig_count):
                f.write('<center><img src=\"section.'+str(i)+ "/fig_"+str(j) + '.png\" width=\"90%\"></center>\n')

            f.write("Comments: [NO COMMENT ADDED]\n<br><br>")

        f.write("</body>\n")
        f.write("</html>\n")
        f.close()




