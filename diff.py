
import phreeqpy.iphreeqc.phreeqc_dll as phreeqc_mod
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tikzplotlib


phreeqc = phreeqc_mod.IPhreeqc("/usr/local/lib/libiphreeqc.so")


def p_initial_solutions(V_w, n_cells, C_0):
    solutions = (
        """
SOLUTION 0
Cl %f
-units 	mol/kgw
-pH 	7	charge_balance
"""
        % C_0
        + "water %.4f# kg\n" % V_w
        + "END\n"
        + """

SOLUTION 1-%d
     -units 	mol/kgw
     -pH 	7	charge_balance\n"""
        % n_cells
        + "water %.4f # kg\n" % V_w
        + "END\n\n"
    )
    return solutions


def p_define_transport(n_cells, steps, dt, dx, d_f):
    transport = (
        "TRANSPORT\n"
        + "cells %d\n" % n_cells
        + "shifts %d\n" % steps
        + f"time_step {dt}\n"
        + f"length {dx}\n"
        + "flow_direction diffusion_only \n"
        + " boundary_conditions constant \n"
        + f"diffusion_coefficient {d_f}\n"
    )
    return transport


def p_selected_output():
    output = """

SELECTED_OUTPUT
    -step    true
    -h    false
    -sim    false
    -soln   true
    -pe     false
    -totals false
    -time   true
    -tot  Cl
    -dist false
    -state false
END
"""
    return output


def p_define_model(n_cells, C_0, d_f, V_w, V_cell, steps, dt, dx):
    model = (
        p_initial_solutions(V_w=V_w, n_cells=n_cells, C_0=C_0)
        + p_define_transport(
            n_cells=n_cells, steps=steps, dt=dt, dx=dx, d_f=d_f
        )
        + p_selected_output()
    )
    return model


n_cells = 20
C_0 = 1
d_f = 5e-10
steps = 3
dt = 1.555e8 / steps
dx = 1 / n_cells
V_cell = 1000 / n_cells  # ml
V_w = 3.8e-1 / n_cells  # kgW

# sigma = 100 * (2 * d_f * dt * steps) ** 0.5
np.set_printoptions(precision=2)
pd.set_option("display.precision", 3)

phreeqc.load_database(
    r"/home/amin/phreeqc_databases/Amm_biozement.dat"
)

phreeqc.run_string(
    p_define_model(
        n_cells=n_cells,
        C_0=C_0,
        d_f=d_f,
        V_w=V_w,
        V_cell=V_cell,
        steps=steps,
        dt=dt,
        dx=dx,
    )
)
p_output = phreeqc.get_selected_output_array()
table_p = pd.DataFrame(p_output[1:], columns=p_output[0])
print("\n")
print(table_p)


row_0 = n_cells + 2  # cell number


def t_plot(var):
    print("row_0= ", row_0)
    label_list = [("k--", 600), ("k-.", 1200), ("k-", 1800)]
    fig, ax = plt.subplots()
    for shift in range(0, steps):
        ax.plot(
            np.arange(0, n_cells) * dx * 100,
            table_p[var][
                row_0
                + shift * n_cells : row_0
                + shift * n_cells
                + n_cells
            ],
            label_list[shift][0],
            label=f"profile after {label_list[shift][1]} days",
            linewidth=2,
        )
        sigma = 100 * (2 * d_f * dt * (shift + 1)) ** 0.5
        ax.plot(
            [2 * sigma, 2 * sigma],
            [0, 1],
            label_list[shift][0],
            label=f"length after {label_list[shift][1]} days",
        )
    ax.legend()
    ax.set(
        xlabel="length (cm)",
        ylabel="C/C$_0$ (relative concentration)",
        title="Diffusion profile and diffusion length",
    )
    ax.set_facecolor("green")
    # fig.set_facecolor("grey")
    tikzplotlib.save("diff_plot.tex")
    plt.savefig("Cl_plot.png")


t_plot(var="Cl(mol/kgw)")


def sigma_t_plot():
    fig, ax = plt.subplots()
    t = np.arange(0, dt * steps, 1e5)
    ax.plot(t / 24 / 3600, 100 * 2 * (2 * d_f * t) ** 0.5, "k")
    ax.set(xlabel="days", ylabel="diffusion length (cm)")
    tikzplotlib.save("sigma_t_plot.tex")


sigma_t_plot()
