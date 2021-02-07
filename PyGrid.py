import argparse
import output as out
from PyGrid_exec import execPyGrid

parser = argparse.ArgumentParser()
parser.add_argument("-Nx", "--Nx", type=int, help="number of elements in x-direction")
parser.add_argument("-Ny", "--Ny", type=int, help="number of elements in y-direction")
parser.add_argument("-fbody","--fbody", type=int, help="geometry type")
parser.add_argument("-plot","--PLOT", action='store_true', help="plot")
parser.add_argument("-q","--Q", type = int, default = 1, help="q is the degree of the curved boundary")
args = parser.parse_args()
Nx = args.Nx
Ny = args.Ny
plot_flag = args.PLOT
q = args.Q
bid = args.fbody

vertices, cell_list, domain = execPyGrid(Nx, Ny, plot_flag, q, bid)
out.output_ply(vertices, cell_list, domain, Nx, Ny, q)
