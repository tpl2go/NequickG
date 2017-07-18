from ComputeModip import *
import matplotlib.pyplot as plt

def visualise_modip_file():
    data = read_stMoDIP('/home/tpl/Documents/Airbus/Project/Papers/Nequick/CCIR_MoDIP/modipNeQG_wrapped.txt')
    plt.contour(data)
    plt.show()

visualise_modip_file()
