import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np


def plot_moments(sols):

    fig, ax = plt.subplots(1, sols[0].state.Nmom)
    for sol in sols:
        for i in range(sol.state.Nmom):
            ax[i].plot(
                sol.times,
                sol.moms[:, i],
                label="NR0 = " + str(sol.state.NR0) + " Nfilt = " + str(sol.Nfilt),
            )
            ax[i].set(xlabel="$t$", ylabel="$M$" + str(sol.state.moments[i]))
            ax[i].legend(loc="upper right")

    plt.tight_layout()


def plot_integrands(sols):

    a_min = 0
    a_max = len(sols[0].times) - 1
    a_init = int(a_max / 2.0)

    fig = plt.figure(figsize=(8, 3))
    plt_ax = plt.axes([0.1, 0.2, 0.8, 0.65])
    slider_ax = plt.axes([0.1, 0.05, 0.8, 0.05])

    plt.axes(plt_ax)

    lin = False
    if abs(sols[0].wave.amplitude - 1) < 1e-2:
        lin = True

    y_dat = sols[0].save[a_init, :, 0]
    if lin:
        y_dat /= sols[0].state.R0

    (my_plot,) = plt.plot(sols[0].state.R0, y_dat)

    a_slider = Slider(slider_ax, "$t$", a_min, a_max, valinit=a_init)

    def update(a):
        y_dat = sols[0].save[int(a), :, 0]
        if lin:
            y_dat /= sols[0].state.R0
            my_plot.set_ydata(y_dat)
        fig.canvas.draw_idle()

    a_slider.on_changed(update)
