"""Transport of a single component in a column with advection and dispersion."""

import phreeqpy.iphreeqc.phreeqc_dll as phreeqc_mod
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tikzplotlib


phreeqc = phreeqc_mod.IPhreeqc("/usr/local/lib/libiphreeqc.so")


def initial_solutions(V_w, n_cells, C_0):
    solutions = (
        f"""
SOLUTION 0
Cl {C_0}
-units 	mol/kgw
-pH 	7	charge_balance\n
"""
        + f"water {V_w}# kg\n"
        + "END\n"
        + """

SOLUTION 1-%d
     -units 	mol/kgw
     -pH 	7	charge_balance\n"""
        % n_cells
        + f"water {V_w} # kg\n"
        + "END\n\n"
    )
    return solutions


def define_transport(n_cells, steps, dt, dx, eps_0, d_f, alpha):
    transport = (
        "TRANSPORT\n"
        + "cells %d\n" % n_cells
        + "shifts %d\n" % steps
        + f"time_step {dt}\n"
        + f"length {dx}\n"
        + "porosities 10*%.4f\n" % eps_0
        + f"diffc {d_f}\n"
        + f"dispersivities 10*{alpha}\n"
        + "flow_direction forward \nboundary_conditions constant \n"
    )
    return transport


def selected_output():
    output = """

SELECTED_OUTPUT
    -step    true
    -h    false
    -sim    false
    -soln   true
    -pe     false
    -totals false
    -time   true
    -pH true
    -tot  Cl
    -dist false
    -state false
END
"""
    return output


def define_model(
    n_cells,
    C_0,
    eps_0,
    d_f,
    V_w,
    steps,
    dt,
    dx,
    alpha,
):
    model = (
        initial_solutions(V_w=V_w, n_cells=n_cells, C_0=C_0)
        + define_transport(
            n_cells=n_cells,
            steps=steps,
            dt=dt,
            dx=dx,
            eps_0=eps_0,
            d_f=d_f,
            alpha=alpha,
        )
        + selected_output()
    )
    return model


# input parameters
n_cells = 40
C_0 = 1
eps_0 = 0.3
d_f = 2e-9  # diffusion coefficient
steps = n_cells  # int(n_cells*1.1)
dt = 1e7 / n_cells
dx = 20 / n_cells
alpha = 2e-1  # dispersivity
V_w = 1  # kgw
u_w = dx / dt  # velocity of water
peclet = u_w * alpha / d_f

print("u_w: ", u_w, " m/s")
np.set_printoptions(precision=2)
pd.set_option("display.precision", 3)

phreeqc.load_database(
    r"/home/amin/phreeqc_databases/Amm_biozement.dat"
)

phreeqc.run_string(
    define_model(
        n_cells=n_cells,
        C_0=C_0,
        eps_0=eps_0,
        d_f=d_f,
        V_w=V_w,
        steps=steps,
        dt=dt,
        dx=dx,
        alpha=alpha,
    )
)
output = phreeqc.get_selected_output_array()
table = pd.DataFrame(output[1:], columns=output[0])
print("\n")
print(table)


row_0 = n_cells + 2  # cell number


def t_plot(var, shift):
    print("row_0= ", row_0)
    fig, ax = plt.subplots()
    ax.plot(
        np.arange(1, n_cells + 1) * dx * 100,
        table[var][
            row_0 + shift * n_cells : row_0 + n_cells * (shift + 1)
        ],
        "ko",
        label=var,
        linewidth=2,
    )
    ax.legend()
    ax.set_facecolor("green")
    fig.set_facecolor("grey")
    plt.savefig("Cl_plot.png")


# t_plot(var="Cl(mol/kgw)", shift=4 * int(n_cells / 10))


def t_plot_mul(var):
    print("row_0= ", row_0)
    fig, ax = plt.subplots()
    label_color = ["k", "b", "red", "orange", "teal"]
    i = 0
    for shift in range(0, steps, int(steps / 5)):
        ax.plot(
            np.arange(1, n_cells + 1) * dx,
            table[var][
                row_0
                + shift * n_cells : row_0
                + n_cells * (shift + 1)
            ],
            label_color[i],
            label=str(shift / n_cells) + " pore volume",
            # + f"{d_f + alpha * u_w:.2e}",
            linewidth=2,
        )
        i += 1
    ax.legend()
    ax.set(
        title=f" Peclet number: {peclet:.2f}, velocity: {u_w:.1e} m/s ",
        xlabel="length (m)",
        ylabel="relative concentration",
    )
    # ax.set_facecolor("grey")
    # fig.set_facecolor("grey")
    plt.savefig("Cl_plot.png")
    tikzplotlib.save("high_peclet_2.tex")


t_plot_mul(var="Cl(mol/kgw)")

# plt.show()
