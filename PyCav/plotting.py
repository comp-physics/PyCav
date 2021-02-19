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
                label="NR0 = " + str(sol.state.NR0) + " Nfilt = " + str(sol.state.Nfilt),
            )
            ax[i].set(xlabel="$t$", ylabel="$M$" + str(sol.state.moments[i]))
            ax[i].legend(loc="lower center",bbox_to_anchor=(0.5,1.03))

    plt.tight_layout()
    plt.show()


def plot_integrands_sub(sols):

    t_min = 0
    t_max = len(sols[0].times) - 1
    t_init = int(t_max / 2.0)

    # fig = plt.figure(figsize=(8, 8))
    fig, ax = plt.subplots(nrows=2, ncols=1)
    plt.subplots_adjust(bottom=0.25)

    slider_ax = plt.axes([0.1, 0.05, 0.8, 0.05])

    lin = is_lin(sols[0].wave.amplitude)
    line0 = [0] * len(sols)
    line1 = [0] * len(sols)
    for i, sol in enumerate(sols):
        x_dat = sol.state.R0
        y_dat = get_RV(t_init,sol,lin)
        # for i in range(2):
        #     ax[i].plot(x_dat,y_dat[i])
        line0[i], = ax[0].plot(x_dat,y_dat[0])
        line1[i], = ax[1].plot(x_dat,y_dat[1])
        for a in ax:
            a.set_xscale('log')

    ax[0].set(ylabel="$f(R_o) R(R_o)$")
    ax[1].set(xlabel="$R_o$", ylabel="$f(R_o) \dot{R}(R_o)$")

    a_slider = Slider(slider_ax, "$t$", t_min, t_max, valinit=t_init)

    def update(t):
        for i, sol in enumerate(sols):
            x_dat = sol.state.R0
            y_dat = get_RV(t,sol,lin)
            line0[i].set_ydata(y_dat[0])
            line1[i].set_ydata(y_dat[1])
        fig.canvas.draw_idle()
 
    a_slider.on_changed(update)
    plt.show()

def plot_integrands(sols):

    t_min = 0
    t_max = len(sols[0].times) - 1
    t_init = int(t_max / 2.0)

    fig = plt.figure(figsize=(8, 3))
    plt_ax = plt.axes([0.1, 0.2, 0.8, 0.65])
    slider_ax = plt.axes([0.1, 0.05, 0.8, 0.05])

    plt.axes(plt_ax)

    lin = is_lin(sols[0].wave.amplitude)
    for sol in sols:
        x_dat = sol.state.R0
        y_dat = get_y(t_init,sol,lin)
        (my_plot,) = plt.plot(x_dat, y_dat)
    plt.xscale("log")

    a_slider = Slider(slider_ax, "$t$", t_min, t_max, valinit=t_init)

    def update(t):
        plt.cla()
        for sol in sols:
            x_dat = sol.state.R0
            y_dat = get_y(t,sol,lin)
            (my_plot,) = plt.plot(x_dat, y_dat)
        plt.xscale("log")
        fig.canvas.draw_idle()
 
    a_slider.on_changed(update)
    plt.show()

def is_lin(p):
    if abs(p - 1) < 1e-2:
        lin = True
    else:
        lin = False
    return lin

def get_y(t=0,sol=None,lin=True):
    p = sol.wave.amplitude
    eps = -1*(p-1.)
    pdf = sol.state.f.copy()
    # pdf = 1.

    y_dat = sol.save[int(t), :, 0]
    if lin:
        y_dat /= sol.state.R0
        y_dat -= 1.
        y_dat /= eps
    y_dat *= pdf
    return y_dat

def get_RV(t=0,sol=None,lin=True):
    p = sol.wave.amplitude
    eps = -1*(p-1.)
    pdf = sol.state.f.copy()
    # pdf = 1.

    R_dat = sol.save[int(t), :, 0]
    V_dat = sol.save[int(t), :, 1]
    if lin:
        R_dat /= sol.state.R0
        R_dat -= 1.
        R_dat /= eps
    R_dat *= pdf
    V_dat *= pdf

    return [R_dat, V_dat]
