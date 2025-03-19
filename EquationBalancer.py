import tkinter as tk
from tkinter import messagebox
import re
import sympy
from math import gcd
import cairo
import io
from PIL import Image, ImageTk

def parse_formula(formula):
    pattern = r'(\([A-Za-z0-9]+\)\d*)|([A-Z][a-z]?\d*)'
    tokens = re.findall(pattern, formula)
    pieces = [t[0] if t[0] else t[1] for t in tokens]
    element_dict = {}
    for piece in pieces:
        if piece.startswith("("):
            inner, subscript = re.findall(r'\(([A-Za-z0-9]+)\)(\d*)', piece)[0]
            subscript = int(subscript) if subscript else 1
            inner_counts = parse_formula(inner)
            for elem, cnt in inner_counts.items():
                element_dict[elem] = element_dict.get(elem, 0) + cnt * subscript
        else:
            match = re.match(r'([A-Z][a-z]?)(\d*)', piece)
            if match:
                elem, num = match.groups()
                num = int(num) if num else 1
                element_dict[elem] = element_dict.get(elem, 0) + num
    return element_dict

def parse_side(side_text):
    compounds = [s.strip() for s in side_text.split('+')]
    return [parse_formula(c) for c in compounds if c]

def parse_equation(left_text, right_text):
    left_compounds = parse_side(left_text)
    right_compounds = parse_side(right_text)
    return left_compounds, right_compounds

def build_matrix(left_compounds, right_compounds):
    all_elements = set()
    for c in left_compounds + right_compounds:
        all_elements.update(c.keys())
    all_elements = sorted(all_elements)
    rows = []
    for element in all_elements:
        row = []
        for comp in left_compounds:
            row.append(comp.get(element, 0))
        for comp in right_compounds:
            row.append(-comp.get(element, 0))
        rows.append(row)
    return sympy.Matrix(rows), all_elements

def solve_coefficients(left_compounds, right_compounds):
    mat, elements = build_matrix(left_compounds, right_compounds)
    ns = mat.nullspace()
    if not ns:
        return None
    sol = ns[0]
    denominators = [x.as_numer_denom()[1] for x in sol]
    lcm_denom = 1
    for d in denominators:
        lcm_denom = sympy.lcm(lcm_denom, d)
    int_sol = [int(lcm_denom * x) for x in sol]
    if any(x < 0 for x in int_sol):
        if all(x <= 0 for x in int_sol):
            int_sol = [-x for x in int_sol]
    g = abs(int_sol[0])
    for v in int_sol[1:]:
        g = gcd(g, abs(v))
    if g != 0:
        int_sol = [v // g for v in int_sol]
    return int_sol

def draw_arrow_line(ctx, x1, y1, x2, y2, arrow_size=8, line_width=2):
    import math
    dx = x2 - x1
    dy = y2 - y1
    angle = math.atan2(dy, dx)
    ctx.save()
    ctx.set_line_width(line_width)
    ctx.set_source_rgb(1, 1, 1)
    ctx.move_to(x1, y1)
    ctx.line_to(x2, y2)
    ctx.stroke()
    ctx.translate(x2, y2)
    ctx.rotate(angle)
    ctx.move_to(0, 0)
    ctx.line_to(-arrow_size, arrow_size / 2)
    ctx.line_to(-arrow_size, -arrow_size / 2)
    ctx.close_path()
    ctx.fill()
    ctx.restore()

def draw_table_cairo_with_arrows(left_compounds, right_compounds, coeffs):
    all_elements = set()
    for c in left_compounds + right_compounds:
        all_elements.update(c.keys())
    all_elements = sorted(all_elements)
    num_left = len(left_compounds)
    num_right = len(right_compounds)
    total_cols = num_left + num_right
    row_count = len(all_elements) + 1
    row_height = 50
    col_width = 180
    top_margin = 20
    left_margin = 20
    width = left_margin + col_width * (total_cols + 1) + left_margin
    height = top_margin + row_height * row_count + top_margin
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
    ctx = cairo.Context(surface)
    ctx.set_source_rgb(0.117, 0.117, 0.184)
    ctx.rectangle(0, 0, width, height)
    ctx.fill()
    def draw_text_centered(cx, cy, text, size=14, bold=False):
        ctx.save()
        ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD if bold else cairo.FONT_WEIGHT_NORMAL)
        ctx.set_font_size(size)
        xbearing, ybearing, twidth, theight, xadv, yadv = ctx.text_extents(text)
        ctx.move_to(cx - twidth/2, cy + theight/2)
        ctx.set_source_rgb(1, 1, 1)
        ctx.show_text(text)
        ctx.restore()
    ctx.set_line_width(2)
    ctx.set_source_rgb(0.5, 0.5, 0.5)
    for r in range(row_count + 1):
        y_line = top_margin + r * row_height
        ctx.move_to(left_margin, y_line)
        ctx.line_to(left_margin + col_width * (total_cols + 1), y_line)
        ctx.stroke()
    for c in range(total_cols + 2):
        x_line = left_margin + c * col_width
        ctx.move_to(x_line, top_margin)
        ctx.line_to(x_line, top_margin + row_height * row_count)
        ctx.stroke()
    header_y = top_margin
    draw_text_centered(left_margin + col_width/2, header_y + row_height/2, "Element", size=14, bold=True)
    for i in range(num_left):
        cx = left_margin + col_width * (i + 1) + col_width/2
        cy = header_y + row_height/2
        draw_text_centered(cx, cy, f"Reactant {i+1}", size=14, bold=True)
    for j in range(num_right):
        cx = left_margin + col_width * (num_left + j + 1) + col_width/2
        cy = header_y + row_height/2
        draw_text_centered(cx, cy, f"Product {j+1}", size=14, bold=True)
    for row_idx, element in enumerate(all_elements):
        y1 = top_margin + row_height * (row_idx + 1)
        y_center = y1 + row_height/2
        draw_text_centered(left_margin + col_width/2, y_center, element, size=13, bold=True)
        for i, comp in enumerate(left_compounds):
            x1 = left_margin + col_width * (i + 1)
            x2 = x1 + col_width
            cx1 = x1 + 10
            cx2 = x2 - 10
            cy = y_center
            base_count = comp.get(element, 0)
            final_count = base_count * coeffs[i]
            if final_count == base_count:
                text_val = str(base_count)
                draw_text_centered((x1 + x2)/2, cy, text_val, size=12, bold=False)
            else:
                text_left = str(base_count)
                text_right = str(final_count)
                ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
                ctx.set_font_size(12)
                _, _, lw, lh, _, _ = ctx.text_extents(text_left)
                _, _, rw, rh, _, _ = ctx.text_extents(text_right)
                ctx.move_to(cx1, cy + lh/2)
                ctx.set_source_rgb(1,1,1)
                ctx.show_text(text_left)
                ctx.move_to(cx2 - rw, cy + rh/2)
                ctx.show_text(text_right)
                arrow_start_x = cx1 + lw + 5
                arrow_end_x = (cx2 - rw) - 5
                arrow_y = cy
                if arrow_end_x > arrow_start_x:
                    draw_arrow_line(ctx, arrow_start_x, arrow_y, arrow_end_x, arrow_y, arrow_size=8, line_width=2)
        for j, comp in enumerate(right_compounds):
            col_index = num_left + j + 1
            x1 = left_margin + col_width * (col_index)
            x2 = x1 + col_width
            cx1 = x1 + 10
            cx2 = x2 - 10
            cy = y_center
            base_count = comp.get(element, 0)
            final_count = base_count * coeffs[num_left + j]
            if final_count == base_count:
                text_val = str(base_count)
                draw_text_centered((x1 + x2)/2, cy, text_val, size=12, bold=False)
            else:
                text_left = str(base_count)
                text_right = str(final_count)
                ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
                ctx.set_font_size(12)
                _, _, lw, lh, _, _ = ctx.text_extents(text_left)
                _, _, rw, rh, _, _ = ctx.text_extents(text_right)
                ctx.move_to(cx1, cy + lh/2)
                ctx.set_source_rgb(1,1,1)
                ctx.show_text(text_left)
                ctx.move_to(cx2 - rw, cy + rh/2)
                ctx.show_text(text_right)
                arrow_start_x = cx1 + lw + 5
                arrow_end_x = (cx2 - rw) - 5
                arrow_y = cy
                if arrow_end_x > arrow_start_x:
                    draw_arrow_line(ctx, arrow_start_x, arrow_y, arrow_end_x, arrow_y, arrow_size=8, line_width=2)
    buf = surface.get_data()
    pil_image = Image.frombuffer("RGBA", (width, height), bytes(buf), "raw", "BGRA", 0, 1)
    return pil_image

class EquationBalancerApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Chemical Equation Balancer")
        self.geometry("1200x600")
        self.configure(bg="#1E1E2F")
        self.create_top_bar()
        self.left_area = tk.Frame(self, bg="#1E1E2F")
        self.left_area.pack(side="left", fill="both", expand=True)
        self.right_sidebar = tk.Frame(self, bg="#222", width=400)
        self.right_sidebar.pack(side="right", fill="y")
        self.input_frame = tk.Frame(self.right_sidebar, bg="#333", height=200)
        self.input_frame.pack(side="top", fill="x")
        self.table_label = tk.Label(self.left_area, bg="#1E1E2F")
        self.table_label.pack(side="top", pady=20)
        self.balanced_equation_label = tk.Label(self.left_area, text="(No results yet)", bg="#1E1E2F", fg="white", font=("Helvetica", 12))
        self.balanced_equation_label.pack(side="top", pady=10)
        self.left_text_var = tk.StringVar()
        self.right_text_var = tk.StringVar()
        self.create_input_section()
        self.table_img_ref = None

    def create_top_bar(self):
        top_bar = tk.Frame(self, bg="#2F416D", height=60)
        top_bar.pack(side="top", fill="x")
        title_label = tk.Label(top_bar, text="Chemical Equation Balancer", fg="white", bg="#2F416D", font=("Helvetica", 18, "bold"))
        title_label.pack(side="left", padx=20)

    def create_input_section(self):
        label_title = tk.Label(self.input_frame, text="Enter Equation", bg="#333", fg="white", font=("Helvetica", 14, "bold"))
        label_title.pack(pady=10)
        reactants_label = tk.Label(self.input_frame, text="Reactants:", bg="#333", fg="white", font=("Helvetica", 12))
        reactants_label.pack(anchor="w", padx=10)
        reactants_entry = tk.Entry(self.input_frame, textvariable=self.left_text_var, font=("Helvetica", 12), width=30)
        reactants_entry.pack(padx=10, pady=5)
        products_label = tk.Label(self.input_frame, text="Products:", bg="#333", fg="white", font=("Helvetica", 12))
        products_label.pack(anchor="w", padx=10)
        products_entry = tk.Entry(self.input_frame, textvariable=self.right_text_var, font=("Helvetica", 12), width=30)
        products_entry.pack(padx=10, pady=5)
        balance_button = tk.Button(self.input_frame, text="Balance", font=("Helvetica", 12, "bold"), bg="#3D5AFE", fg="white", command=self.on_balance_clicked)
        balance_button.pack(pady=10)

    def on_balance_clicked(self):
        left_text = self.left_text_var.get().strip()
        right_text = self.right_text_var.get().strip()
        if not left_text or not right_text:
            messagebox.showwarning("Input Error", "Please enter both reactants and products.")
            return
        try:
            left_compounds, right_compounds = parse_equation(left_text, right_text)
        except Exception as e:
            messagebox.showerror("Parsing Error", f"Could not parse equation:\n{e}")
            return
        if not left_compounds or not right_compounds:
            messagebox.showwarning("Input Error", "Invalid compounds on one or both sides.")
            return
        coeffs = solve_coefficients(left_compounds, right_compounds)
        if coeffs is None:
            messagebox.showerror("Balancing Error", "No non-trivial solution found.")
            return
        num_left = len(left_compounds)
        left_coeffs = coeffs[:num_left]
        right_coeffs = coeffs[num_left:]
        left_side_raw = [s.strip() for s in left_text.split('+')]
        right_side_raw = [s.strip() for s in right_text.split('+')]
        balanced_left = []
        for c, formula_str in zip(left_coeffs, left_side_raw):
            balanced_left.append(formula_str if c == 1 else f"{c}{formula_str}")
        balanced_right = []
        for c, formula_str in zip(right_coeffs, right_side_raw):
            balanced_right.append(formula_str if c == 1 else f"{c}{formula_str}")
        balanced_eq_str = " + ".join(balanced_left) + " → " + " + ".join(balanced_right)
        self.balanced_equation_label.config(text=f"Balanced Equation:\n{balanced_eq_str}")
        pil_image = draw_table_cairo_with_arrows(left_compounds, right_compounds, coeffs)
        self.table_img_ref = ImageTk.PhotoImage(pil_image)
        self.table_label.config(image=self.table_img_ref)

if __name__ == "__main__":
    app = EquationBalancerApp()
    app.mainloop()
