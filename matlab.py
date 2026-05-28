import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt

# =========================
# HELPERS
# =========================

def fmt(v):
    return f"{v:.3f}" if isinstance(v, float) else v

def r(v):
    return round(v, 3)

def f(x, eq):
    return eval(eq, {"x": x, "np": np})

def df(x, eq):
    h = 1e-5
    return (f(x + h, eq) - f(x - h, eq)) / (2 * h)

def stop(ea, val, tol):
    return abs(val) < tol or (ea != 0 and ea < tol)

# =========================
# ROOT FINDING METHODS
# =========================

def bisection(eq, xl, xu, tol, max_iter):

    if xl > xu:
        xl, xu = xu, xl

    if f(xl, eq) * f(xu, eq) > 0:
        raise ValueError("f(XL) and f(XU) must have opposite signs")

    t = []
    xr_old = None

    for i in range(1, max_iter + 1):

        xr = (xl + xu) / 2

        fxl = f(xl, eq)
        fxr = f(xr, eq)

        ea = 0 if xr_old is None else abs((xr - xr_old) / xr) * 100

        prod = fxl * fxr

        rem = "<0 Revert back to XL" if prod < 0 else ">0 Go to next interval"

        t.append((
            i,
            r(xl),
            r(xr),
            r(xu),
            r(fxl),
            r(fxr),
            r(ea),
            r(prod),
            rem
        ))

        if stop(ea, fxr, tol):
            break

        if prod < 0:
            xu = xr
        else:
            xl = xr

        xr_old = xr

    return xr, t


def false_position(eq, xl, xu, tol, max_iter):

    if xl > xu:
        xl, xu = xu, xl

    if f(xl, eq) * f(xu, eq) > 0:
        raise ValueError("f(XL) and f(XU) must have opposite signs")

    t = []
    xr_old = None

    for i in range(1, max_iter + 1):

        fxl = f(xl, eq)
        fxu = f(xu, eq)

        xr = xu - (fxu * (xl - xu)) / (fxl - fxu)

        fxr = f(xr, eq)

        ea = 0 if xr_old is None else abs((xr - xr_old) / xr) * 100

        prod = fxl * fxr

        rem = "<0 Revert back to XL" if prod < 0 else ">0 Go to next interval"

        t.append((
            i,
            r(xl),
            r(xu),
            r(xr),
            r(ea),
            r(fxl),
            r(fxu),
            r(fxr),
            r(prod),
            rem
        ))

        if stop(ea, fxr, tol):
            break

        if prod < 0:
            xu = xr
        else:
            xl = xr

        xr_old = xr

    return xr, t


def secant(eq, x0, x1, tol, max_iter):

    t = []

    for i in range(1, max_iter + 1):

        f0 = f(x0, eq)
        f1 = f(x1, eq)

        if f0 == f1:
            break

        x2 = x1 - f1 * (x0 - x1) / (f0 - f1)

        f2 = f(x2, eq)

        ea = abs((x2 - x1) / x2) * 100

        t.append((
            i,
            r(x0),
            r(x1),
            r(x2),
            r(ea),
            r(f0),
            r(f1),
            r(f2)
        ))

        if stop(ea, f2, tol):
            break

        x0 = x1
        x1 = x2

    return x2, t


def newton(eq, x0, tol, max_iter):

    t = []

    for i in range(1, max_iter + 1):

        fx = f(x0, eq)

        dfx = df(x0, eq)

        if dfx == 0:
            break

        x1 = x0 - fx / dfx

        ea = abs((x1 - x0) / x1) * 100

        t.append((
            i,
            r(x0),
            r(ea),
            r(fx),
            r(dfx)
        ))

        if stop(ea, fx, tol):
            break

        x0 = x1

    return x1, t


def incremental(eq, x0, ΔX, tol, max_iter):

    t = []

    for i in range(1, max_iter + 1):

        x1 = x0 + ΔX

        f0 = f(x0, eq)
        f1 = f(x1, eq)

        prod = f0 * f1

        rem = "<0 Revert back to XL" if prod < 0 else ">0 Go to next interval"

        t.append((
            i,
            r(x0),
            r(ΔX),
            r(x1),
            r(f0),
            r(f1),
            r(prod),
            rem
        ))

        if abs(f1) < tol or prod < 0:
            return (x0, x1), t

        x0 = x1

    return None, t

# =========================
# GRAPH
# =========================

def graph(eq):

    try:

        m = method.get()

        tol      = get_tol()
        max_iter = get_max_iter()

        points = []

        if m in ["Bisection", "False Position"]:

            xl = float(entries["XL"].get())
            xu = float(entries["XU"].get())

            if xl > xu:
                xl, xu = xu, xl

            xmin = xl - 2
            xmax = xu + 2

            if m == "Bisection":
                rootv, _ = bisection(eq, xl, xu, tol, max_iter)
            else:
                rootv, _ = false_position(eq, xl, xu, tol, max_iter)

            xr = rootv

            points.append(("XL", xl, f(xl, eq), "blue"))
            points.append(("XU", xu, f(xu, eq), "green"))
            points.append(("XR", xr, f(xr, eq), "red"))

        elif m == "Secant":

            x0 = float(entries["Xi-1"].get())
            x1 = float(entries["Xi"].get())

            xmin = min(x0, x1) - 2
            xmax = max(x0, x1) + 2

            rootv, _ = secant(eq, x0, x1, tol, max_iter)

            points.append(("Xi-1", x0, f(x0, eq), "blue"))
            points.append(("Xi", x1, f(x1, eq), "green"))
            points.append(("Root", rootv, f(rootv, eq), "red"))

        elif m == "Newton":

            x0 = float(entries["Xi"].get())

            xmin = x0 - 5
            xmax = x0 + 5

            rootv, _ = newton(eq, x0, tol, max_iter)

            points.append(("Xi", x0, f(x0, eq), "blue"))
            points.append(("Root", rootv, f(rootv, eq), "red"))

        else:

            xl = float(entries["XL"].get())
            ΔX = float(entries["ΔX"].get())

            xmin = xl - 2
            xmax = xl + ΔX * 10

            rr, _ = incremental(eq, xl, ΔX, tol, max_iter)

            if rr:
                x0, x1 = rr

                rootv = (x0 + x1) / 2

                points.append(("XL", x0, f(x0, eq), "blue"))
                points.append(("XU", x1, f(x1, eq), "green"))
                points.append(("Root", rootv, f(rootv, eq), "red"))

        x = np.linspace(xmin, xmax, 1000)

        y = [f(i, eq) for i in x]

        plt.figure(figsize=(10, 5))

        plt.plot(x, y, linewidth=2, label="f(x)")

        plt.axhline(0, color='black')

        for name, px, py, color in points:

            plt.scatter(px, py, s=120, color=color, zorder=5)

            plt.text(
                px,
                py,
                f"{name}\n({px:.3f}, {py:.3f})",
                fontsize=9
            )

        plt.grid(True)

        plt.title(f"{m} Method Graph")

        plt.xlabel("x")
        plt.ylabel("f(x)")

        plt.legend()

        plt.show()

    except Exception as e:
        messagebox.showerror("Graph Error", str(e))

# =========================
# MATRIX
# =========================

def matrix(txt):

    try:
        return np.array([
            list(map(float, row.split()))
            for row in txt.strip().split(";")
        ])

    except:
        messagebox.showerror(
            "Error",
            "Format:\n1 2;3 4"
        )

# =========================
# GUI
# =========================

root = tk.Tk()

root.title("Numerical Methods PRO")

root.geometry("1150x720")

tab = ttk.Notebook(root)

tab.pack(fill="both", expand=True)

# =========================
# ROOT FINDING TAB
# =========================

f1 = tk.Frame(tab)

tab.add(f1, text="Root Finding")

tk.Label(f1, text="Equation f(x)").pack()

eq = tk.Entry(f1, width=50)

eq.insert(0, "x**3-x-2")

eq.pack()

method = ttk.Combobox(
    f1,
    values=[
        "Incremental",
        "Bisection",
        "False Position",
        "Newton",
        "Secant"
    ]
)

method.pack()

method.set("Bisection")

# --- Tolerance and Max Iterations row
tol_frame = tk.Frame(f1)
tol_frame.pack()

tk.Label(tol_frame, text="Tolerance:").pack(side="left", padx=(0, 4))
tol_entry = tk.Entry(tol_frame, width=10)
tol_entry.insert(0, "0.001")
tol_entry.pack(side="left", padx=(0, 20))

tk.Label(tol_frame, text="Max Iterations:").pack(side="left", padx=(0, 4))
max_iter_entry = tk.Entry(tol_frame, width=6)
max_iter_entry.insert(0, "50")
max_iter_entry.pack(side="left")

def get_tol():
    try:
        v = float(tol_entry.get())
        return v if v > 0 else 0.001
    except:
        return 0.001

def get_max_iter():
    try:
        v = int(max_iter_entry.get())
        return v if v >= 1 else 50
    except:
        return 50

header = tk.Label(
    f1,
    fg="blue",
    font=("Arial", 10, "bold")
)

header.pack()

inputs = tk.Frame(f1)

inputs.pack()

entries = {}

def update(e=None):

    for w in inputs.winfo_children():
        w.destroy()

    entries.clear()

    m = method.get()

    data = {

        "Bisection": (
            "Iteration | XL | XR | XU | f(XL) | f(XR) | Ea | Product | Remarks",
            ["XL", "XU"]
        ),

        "False Position": (
            "Iteration | XL | XU | XR | Ea | f(XL) | f(XU) | f(XR) | Product",
            ["XL", "XU"]
        ),

        "Secant": (
            "Iteration | Xi-1 | Xi | Xi+1 | Ea",
            ["Xi-1", "Xi"]
        ),

        "Newton": (
            "Iteration | Xi | Ea | f(x) | f'(x)",
            ["Xi"]
        ),

        "Incremental": (
            "Iteration | XL | ΔX | XU | f(XL) | f(XU) | Product | Remarks",
            ["XL", "ΔX"]
        ),
    }

    header.config(text=data[m][0])

    for n in data[m][1]:

        tk.Label(inputs, text=n).pack()

        e = tk.Entry(inputs)

        e.pack()

        entries[n] = e

method.bind("<<ComboboxSelected>>", update)

update()

tree = ttk.Treeview(f1)

tree.pack(fill="both", expand=True)

def solve():

    tree.delete(*tree.get_children())

    try:

        m        = method.get()
        tol      = get_tol()
        max_iter = get_max_iter()

        if m == "Bisection":

            rootv, t = bisection(
                eq.get(),
                float(entries["XL"].get()),
                float(entries["XU"].get()),
                tol,
                max_iter
            )

            cols = (
                "i","XL","XR","XU",
                "f(XL)","f(XR)",
                "Ea","Product","Remarks"
            )

        elif m == "False Position":

            rootv, t = false_position(
                eq.get(),
                float(entries["XL"].get()),
                float(entries["XU"].get()),
                tol,
                max_iter
            )

            cols = (
                "i","XL","XU","XR",
                "Ea","f(XL)",
                "f(XU)","f(XR)",
                "Product","Remarks"
            )

        elif m == "Secant":

            rootv, t = secant(
                eq.get(),
                float(entries["Xi-1"].get()),
                float(entries["Xi"].get()),
                tol,
                max_iter
            )

            cols = (
                "i","Xi-1","Xi",
                "Xi+1","Ea",
                "f(Xi-1)","f(Xi)",
                "f(Xi+1)"
            )

        elif m == "Newton":

            rootv, t = newton(
                eq.get(),
                float(entries["Xi"].get()),
                tol,
                max_iter
            )

            cols = (
                "i","Xi","Ea",
                "f(x)","f'(x)"
            )

        else:  # Incremental

            rootv, t = incremental(
                eq.get(),
                float(entries["XL"].get()),
                float(entries["ΔX"].get()),
                tol,
                max_iter
            )

            cols = (
                "i","XL","ΔX","XU",
                "f(XL)","f(XU)",
                "Product","Remarks"
            )

        tree["columns"] = cols

        tree["show"] = "headings"

        for c in cols:

            tree.heading(c, text=c)

            tree.column(c, width=120)

        for row in t:

            tree.insert(
                "",
                tk.END,
                values=tuple(fmt(v) for v in row)
            )

        # ── FIX: Incremental returns a tuple (x0, x1), handle it separately ──
        if rootv is None:
            messagebox.showinfo("Result", "No root found in the given interval.")
        elif isinstance(rootv, tuple):
            x0, x1 = rootv
            messagebox.showinfo(
                "Result",
                f"Root bracket found!\nXL ≈ {x0:.3f}  |  XU ≈ {x1:.3f}\nEstimated root ≈ {(x0 + x1) / 2:.3f}"
            )
        else:
            messagebox.showinfo(
                "Result",
                f"Root ≈ {rootv:.3f}"
            )

    except Exception as e:

        messagebox.showerror(
            "Error",
            str(e)
        )

tk.Button(
    f1,
    text="Solve",
    bg="green",
    fg="white",
    command=solve
).pack(pady=5)

tk.Button(
    f1,
    text="Plot Graph",
    bg="blue",
    fg="white",
    command=lambda: graph(eq.get())
).pack(pady=5)

# =========================
# MATRIX TAB
# =========================

f2 = tk.Frame(tab)

tab.add(f2, text="Matrix")

tk.Label(
    f2,
    text="Matrix A (example: 1 2;3 4)"
).pack()

A = tk.Entry(f2, width=50)

A.pack()

tk.Label(f2, text="Matrix B").pack()

B = tk.Entry(f2, width=50)

B.pack()

box = tk.Text(f2, height=15)

box.pack()

def show(x):

    box.delete("1.0", tk.END)

    box.insert(tk.END, str(x))

ops = {

    "Add":
        lambda a,b: a+b,

    "Multiply":
        lambda a,b: a@b,

    "Transpose":
        lambda a,b: a.T,

    "Determinant":
        lambda a,b: np.linalg.det(a),

    "Inverse":
        lambda a,b: np.linalg.inv(a),

    "Adjoint":
        lambda a,b: np.linalg.inv(a)*np.linalg.det(a),

    "Power":
        lambda a,b: np.linalg.matrix_power(a,2),

    "Solve Ax=b":
        lambda a,b: np.linalg.solve(
            a,
            np.array([5,6])
        )
}

for n, func in ops.items():

    tk.Button(
        f2,
        text=n,
        width=15,
        command=lambda f=func:
        show(f(matrix(A.get()), matrix(B.get())))
    ).pack()

root.mainloop()