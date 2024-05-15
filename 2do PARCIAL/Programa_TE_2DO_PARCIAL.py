from tkinter import *
from tkinter import messagebox
import sympy as sp
from math import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import matplotlib.pyplot as plt

# Habilitar uso de texto TeX
#plt.rcParams['text.usetex'] = True
#plt.rcParams['font.family'] = 'serif'
#plt.rcParams['font.serif'] = 'Computer Modern'

#def formatear_complejo(sol):
    #"""Formatea un número complejo para LaTeX reemplazando 'j' por 'i' y redondeando las partes."""
    #re, im = sol.as_real_imag()
    #re, im = re.evalf(), im.evalf()
    #return (sp.latex(re.simplify()), sp.latex(im.simplify()))
def formatear_complejo(sol):
    """Formatea un número complejo para LaTeX reemplazando 'j' por 'i' y manteniendo formato fraccional o raíz cuadrada."""
    re, im = sol.as_real_imag()
    return (sp.latex(re.simplify()), sp.latex(im.simplify()))  # Devuelve una tupla con las partes real e imaginaria en formato LaTeX.


def solucion():
    try:
        a, b, c = (sp.nsimplify(float(entrada.get())) for entrada in entradas)
    except ValueError:
        messagebox.showerror("Error", "Introduce valores numéricos válidos.")
        return

    m = sp.symbols('m')
    ecuacion = sp.Eq(a*m**2 + b*m + c, 0)
    soluciones = sp.solve(ecuacion, m)
    discriminante = b**2 - 4*a*c

# Formatear la ecuación característica
    ec_homo = f"{sp.latex(a)}m^2 + {sp.latex(b)}m + {sp.latex(c)} = 0"
    sol_text = f"Ecuación Característica: ${ec_homo}$\n"

    # Determinar tipo de raíces y formatear soluciones
    if discriminante > 0:
        roots = [formatear_complejo(sol) for sol in soluciones]
        sol_text += "Raíces reales y distintas: " + ", ".join(f"$m = {root[0]}$" for root in roots) + "\n"
        sol_general = " + ".join(f"$C_{i+1}$ $e^{{{root[0]}x}}$" for i, root in enumerate(roots))
    elif discriminante == 0:
        root = formatear_complejo(soluciones[0])
        sol_text += f"Raíces reales y repetidas: $m = {root[0]}$\n"
        sol_general = f"$C_1$ $e^{{{root[0]}x}}$ + $C_2 x$ $e^{{{root[0]}x}}$"
    else:
        roots = [formatear_complejo(sol) for sol in soluciones]
        sol_text += "Raíces complejas: " + ", ".join(f"$m =$ ${root[0]} + {root[1]}$$i$" for root in roots) + "\n"
        sol_general = " + ".join(f"$C_{2*i+1}$$\\cos({root[1]}x)$$e^{{{root[0]}x}}$ +\n $C_{2*i+2}$$\\sin({root[1]}x)$$e^{{{root[0]}x}}$ " for i, root in enumerate(roots))

    sol_text += f"\nSolución General: $Y_C$ =  {sol_general}"

    # El resto del código para la interfaz gráfica se mantiene igual.

    fig = Figure(figsize=(9, 4))
    ax = fig.add_subplot(111)
    ax.axis('off')
    ax.text(0.5, 0.5, sol_text, fontsize=12, va='center', ha='center', wrap=True)

    canvas = FigureCanvasTkAgg(fig, master=miFrame)
    canvas.draw()
    canvas.get_tk_widget().grid(row=4, column=0, columnspan=8)

raiz = Tk()
raiz.title("Solucionador de Ecuaciones Diferenciales")
raiz.iconbitmap("icono.ico")
raiz.config(bg="#B0DDD8")

miFrame = Frame(raiz, width="850", height="550", bg="#ACD3CF")
miFrame.pack_propagate(True)
miFrame.pack()

miLabel = Label(miFrame, text="Solucionador de Ecuaciones Dif. homogéneas de 2do Orden", font=("SF Pro Display", 18), bg="#ACD3CF")
miLabel.grid(column=0, row=0, columnspan=8, padx=10, pady=20)

Label(miFrame, text="Introduce los coeficientes de tu ecuación:", font=("SF Pro Display", 15, "italic"), bg="#ACD3CF").grid(column=0, row=1, columnspan=7, sticky="w", padx=15, pady=5)

entradas = [Entry(miFrame, width=10) for _ in range(3)]
labels = ["y'' +", "y' +", "y = 0"]
for i, entrada in enumerate(entradas):
    entrada.grid(row=3, column=2*i+1, padx=2.5, pady=15)
    Label(miFrame, text=labels[i], fg="#282828", font=("SF Pro Display", 16), bg="#ACD3CF").grid(row=3, column=2*i+2)

boton_calcula = Button(miFrame, text="RESOLVER", command=solucion, font=("SF Pro Display", 13), bg="#65B8B0")
boton_calcula.grid(row=3, column=7, padx=10, pady=10)

raiz.mainloop()
