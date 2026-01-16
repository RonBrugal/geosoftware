import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext
import math
import numpy as np
from typing import List, Dict, Tuple, Optional


def convert_spt_to_n60(n_field, hammer_efficiency=0.6, borehole_diameter=1.0, sampler_correction=1.0,
                       rod_length_correction=1.0):
    """Convert SPT field N-value to N60 using energy and equipment corrections."""
    # Energy ratio correction
    # N60 assumes 60% efficiency
    ER = hammer_efficiency / 0.6

    # Apply all correction factors
    # N60 = N-field × ER × CB × CS × CR
    N60 = n_field * ER * borehole_diameter * sampler_correction * rod_length_correction
    return N60

def convert_spt_to_n160(n_field, hammer_efficiency=0.6, overburden_pressure=0, atmospheric_pressure=2116.0,
                        borehole_diameter=1.0, sampler_correction=1.0, rod_length_correction=1.0):
    """Convert SPT field N-value to N160 (N1)60 using energy, equipment, and overburden correction."""
    # First calculate N60 with all equipment corrections
    N60 = convert_spt_to_n60(n_field, hammer_efficiency, borehole_diameter, sampler_correction, rod_length_correction)

    # Overburden correction factor CN
    # CN = (Pa / σ'v)^0.5 where Pa is atmospheric pressure and σ'v is effective overburden
    CN = (atmospheric_pressure / overburden_pressure) ** 0.5

    N160 = N60 * CN
    return N160

class MatFoundationSettlement:
    """
    Calculate settlement under a mat foundation using Boussinesq's theory
    and geotechnical soil parameters.
    """

    def __init__(self, width, length, depth, borings_data):
        """
        Initialize the mat foundation settlement calculator.

        Parameters:
        -----------
        width : float
            Width of mat foundation (ft)
        length : float
            Length of mat foundation (ft)
        depth : float
            Depth of mat foundation below ground surface (ft)
        borings_data : dict
            Dictionary with boring_id as key and dict containing 'depth', 'E' and 'n60' lists
            Example: {'B-1': {'depth': [5, 10, 15], 'n60': [10, 15, 20]}}
            Example: {'B-1': {'E': [1000, 1500, 2000]}}
        """
        self.width = width
        self.length = length
        self.depth = depth
        self.borings_data = borings_data

        # Storage for loading conditions
        self.point_load = 0  # lbs
        self.uniform_load = 2000  # psf

    def add_point_load(self, p):
        """
        Add a point load to the center of the mat foundation.

        Parameters:
        -----------
        P : float
            Point load magnitude (lbs)
        """
        self.point_load = p

    def set_uniform_load(self, q):
        """
        Set uniform distributed load on the mat foundation.

        Parameters:
        -----------
        q : float
            Uniform distributed load (psf)
        """
        self.uniform_load = q

    def _boussinesq_total_load(self, q, B, L, z):
        """
        Calculate vertical stress at depth z below center of rectangular area
        with uniform load using influence factor approximation.

        Parameters:
        -----------
        q : float
            Uniform load (psf)
        B : float
            Width of loaded area (ft)
        L : float
            Length of loaded area (ft)
        z : float
            Depth below foundation (ft)

        Returns:
        --------
        sigma_z : float
            Vertical stress at depth z (psf)
        """
        if z <= 0:
            return q

        # Simplified influence factor for rectangular loaded area
        # Using approximate formula for center point
        m1 = L / B
        n1 = z / (B / 2)
        #Influence factor I at center of rectangular area (Fadum's chart approximation)

        term1 = (m1 * n1) / (math.sqrt(1 + m1 ** 2 + n1 ** 2) * (1 + n1 ** 2) * (m1 ** 2 + n1 ** 2))
        term2 = (1 + m1 ** 2 + 2 * n1 ** 2) / ((1 + n1 ** 2) * (m1 ** 2 + n1 ** 2))
        term3 = math.asin(m1 / (math.sqrt(m1 ** 2 + n1 ** 2) * math.sqrt(1 + n1 ** 2)))
        I_c = (2 / math.pi) * (term1 * term2 + term3)
        sigma_z_pl = self.point_load / (self.length * self.width)
        sigma_z_tot = (q+sigma_z_pl) * I_c
        return sigma_z_tot

    def calculate_stress_under_footing(self, boring_id):
        """
        Calculate stress at all depths.

        Parameters:
        -----------
        boring_id : str
            Boring identifier

        Returns:
        --------
        stress_profile : list
            List of stresses at each sampled depth (psf)
        """
        if boring_id not in self.borings_data:
            raise ValueError(f"Boring ID '{boring_id}' not found in borings_data")

        depths = self.borings_data[boring_id]['depths']
        stress_profile = []

        for depth in depths:
            # Adjust depth to account for foundation depth
            z = depth - self.depth

            if z <= 0:
                stress_profile.append(0)
                continue

            # Calculate stress from uniform load
            sigma_total = self._boussinesq_total_load(
                self.uniform_load, self.width, self.length, z)

            # Total stress
            stress_profile.append(sigma_total)

        return stress_profile

    def calculate_settlement(self, boring_id):
        """
        Calculate immediate settlement.

        Parameters:
        -----------
        boring_id : str
            Boring identifier
        Returns:
        --------
        settlement_profile : list
            List of tuples (depth, stress, n60, E, settlement) for each layer
        """
        if boring_id not in self.borings_data:
            raise ValueError(f"Boring ID '{boring_id}' not found in borings_data")
        depths = self.borings_data[boring_id]['depths']
        n60_values = self.borings_data[boring_id]['n60_values']
        E_values = self.borings_data[boring_id]['E']
        stress_profile = self.calculate_stress_under_footing(boring_id)

        settlement_profile = []

        for i, (depth, stress, n60, E) in enumerate(zip(depths, stress_profile, n60_values, E_values)):

            # Determine layer thickness (distance to next depth or assumed increment)
            if i < len(depths) - 1:
                layer_thickness = depths[i + 1] - depth
            else:
                layer_thickness = 5  # Assume 5 ft for last layer

            settlement = ((stress * layer_thickness) / (E*1000))/12

            settlement_profile.append({
                'depth': depth,
                'stress': stress,
                'n60': n60,
                'E': E,
                'settlement': settlement
            })

        return settlement_profile

    def get_total_settlement(self, boring_id):
        """
        Get total settlement by summing all layers.

        Parameters:
        -----------
        boring_id : str
            Boring identifier

        Returns:
        --------
        total_settlement : float
            Total immediate settlement (ft)
        """
        settlement_profile = self.calculate_settlement(boring_id)
        total_settlement = sum(layer['settlement'] for layer in settlement_profile)
        return total_settlement

class DrainedModulus:
    def __init__(self, root, borings_data=None):
        self.root = root
        self.root.title("Drained Modulus Stiffness Calculator")
        self.root.geometry("800x700")

        self.borings_data = borings_data
        # Soil type to formula mapping
        # Each soil type has a different formula for calculating E from N60
        self.soil_formulas = {
            1: {
                'name': 'Sand (NC)',
                'formula': 'E = 500(N60 + 15) in kPa',
                'description': 'Bowles 1996'
            },
            2: {
                'name': 'Sand (saturated)',
                'formula': 'E = 250(N60 + 15) in Kpa',
                'description': 'Bowles 1996'
            },
            3: {
                'name': 'Sand (OC)',
                'formula': 'E = 40000 + 1050N60 in kPa',
                'description': 'Bowles 1996'
            },
            4: {
                'name': 'Gravelly Sand',
                'formula': 'E = 600(N60 + 6) for N<= 15\nE = 600(N + 6) + 2000 for N> 15',
                'description': 'Bowles 1996'
            },
            5: {
                'name': 'Clayey Sand',
                'formula': 'E = 320(N60 + 15)',
                'description': 'Bowles 1996'
            },
            6: {
                'name': 'Silts, sandy silts or clayey silt',
                'formula': 'E = 300(N60 + 6',
                'description': 'Bowles 1996'
            }
        }

        # Create main frame
        main_frame = ttk.Frame(root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Configure grid weights
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)

        # Title
        title = ttk.Label(main_frame, text="Drained Modulus Stiffness Calculator", font=('Arial', 16, 'bold'))
        title.grid(row=0, column=0, columnspan=2, pady=10)

        # Help text
        help_text = "Calculates E (kPa) using soil-specific formulas with N60 as input"
        ttk.Label(main_frame, text=help_text, font=('Arial', 9), foreground='gray').grid(row=1, column=0, columnspan=2,
                                                                                         pady=5)

        # Buttons frame
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=2, column=0, columnspan=2, pady=15)
        ttk.Button(button_frame, text="Calculate E", command=self.calculate_modulus).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="View Formulas", command=self.show_formulas).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Back to Main Menu", command=self.back_to_main).pack(side=tk.LEFT, padx=5)

        self.results_text = scrolledtext.ScrolledText(main_frame, width=80, height=20, wrap=tk.WORD)
        self.results_text.grid(row=3, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)

        # Configure text tags for formatting
        self.results_text.tag_configure("header", font=('Courier', 10, 'bold'))
        self.results_text.tag_configure("data", font=('Courier', 10))

        main_frame.rowconfigure(3, weight=1)

    def back_to_main(self):
        self.root.destroy()
        main_menu()

    def calculate_e_from_formula(self, soil_type, n60):
        """
        Calculate E (ksf) based on soil type and N60 value.
        Each soil type uses a different formula.
        """
        if soil_type == 1:  # Sand (NC)
            E = 10*(n60+15)
            return E

        elif soil_type == 2:  # Sand (saturated)
            E = 5*(n60+15)
            return E

        elif soil_type == 3:  # Sand (OC)
            E = (40000 + 1050*n60)/50
            return E

        elif soil_type == 4:  # Gravelly Sand
            if n60 > 15:
                E = (((n60 + 6) * 600 )+ 2000) / 50
            else:
                E = ((n60 + 6) * 600) / 50
            return E

        elif soil_type == 5:  # Clayey Sand
            E = (320*(n60+15))/50
            return E

        elif soil_type == 6:  # Silts, sandy silts or clayey silt
            E = (300*(n60+6)) / 50
            return E

        else:
            return None

    def show_formulas(self):
        """Display the formulas used for each soil type"""
        info_text = "Formulas for Calculating E (kPa) from N60:\n\n"
        info_text += "=" * 70 + "\n"

        for soil_type, data in sorted(self.soil_formulas.items()):
            info_text += f"\nSoil Type {soil_type}: {data['name']}\n"
            info_text += f"Formula: {data['formula']}\n"
            info_text += f"Description: {data['description']}\n"
            info_text += "-" * 70 + "\n"

        messagebox.showinfo("Soil Type Formulas", info_text)

    def calculate_modulus(self):
        if not self.borings_data:
            messagebox.showerror("Error", "Please add boring data with soil types first")
            return

        self.results_text.delete(1.0, tk.END)

        for boring_id, data in self.borings_data.items():
            if 'n60_values' not in data:
                self.results_text.insert(tk.END,
                                         f"Error: No N60 values found for {boring_id}. Run N60 calculation first.\n",
                                         "header")
                continue
            self.results_text.insert(tk.END, f"\n{'=' * 85}\n", "header")
            self.results_text.insert(tk.END, f"BORING: {boring_id}\n", "header")
            self.results_text.insert(tk.END, f"{'=' * 85}\n", "header")

            header = f"{'Depth (ft)':<12} {'Soil Type':<12} {'Soil Name':<25} {'N60':<12} {'E (kPa)':<12}\n"
            self.results_text.insert(tk.END, header, "header")
            self.results_text.insert(tk.END, "-" * 85 + "\n", "header")

            for depth, soil_type, n60 in zip(data['depths'], data['soil types'], data['n60_values']):
                if soil_type in self.soil_formulas:
                    soil_name = self.soil_formulas[soil_type]['name']
                    E = self.calculate_e_from_formula(soil_type, n60) if n60 is not None else None

                    # store E back into borings_data
                    if 'E' not in data:
                        data['E'] = []

                    data['E'].append(E)

                    if E is not None:
                        line = f"{depth:<12.1f} {soil_type:<12} {soil_name:<25} {n60:<12.2f} {E:<12.2f}\n"
                    else:
                        line = f"{depth:<12.1f} {soil_type or 'N/A':<12} {soil_name:<25} {'N/A':<12} {'N/A':<12}\n"

                    self.results_text.insert(tk.END, line, "data")
                else:
                    line = f"{depth:<12.1f} {soil_type:<12} {'Unknown':<25} {n60:<12.2f} {'N/A':<12}\n"
                    self.results_text.insert(tk.END, line, "data")
            self.results_text.insert(tk.END, "\n")

        # Add summary
        self.results_text.insert(tk.END, f"\n{'=' * 85}\n", "header")
        self.results_text.insert(tk.END, f"CALCULATION NOTES:\n", "header")
        self.results_text.insert(tk.END, f"{'=' * 85}\n", "header")
        self.results_text.insert(tk.END, f"Each soil type uses a different formula to calculate E from N60\n", "data")
        self.results_text.insert(tk.END, f"Click 'View Formulas' button to see the formulas for each soil type\n",
                                 "data")
        mat_calcs = MatFoundationSettlement(width=5, length=5, depth=2, borings_data=self.borings_data)
        for boring_id in self.borings_data.keys():
            print(mat_calcs.calculate_settlement(boring_id=boring_id))

class N60:
    def __init__(self, root, borings_data=None):
        self.root = root
        self.root.title("N60 and N160 Calculator")
        self.root.geometry("900x800")

        self.borings_data = borings_data

        # Create main frame
        main_frame = ttk.Frame(root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Configure grid weights
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)

        # Title
        title = ttk.Label(main_frame, text="N60 and N160 Calculator", font=('Arial', 16, 'bold'))
        title.grid(row=0, column=0, columnspan=2, pady=10)

        # Correction parameters section
        params_label = ttk.Label(main_frame, text="Energy Correction", font=('Arial', 12, 'bold'))
        params_label.grid(row=1, column=0, columnspan=2, pady=5)

        # Hammer efficiency
        ttk.Label(main_frame, text="Hammer Efficiency (decimal):").grid(row=2, column=0, sticky=tk.W, pady=5)
        self.hammer_eff_var = tk.StringVar(value="0.80")
        ttk.Entry(main_frame, textvariable=self.hammer_eff_var, width=10).grid(row=2, column=1, sticky=tk.W, pady=5)

        # Help text
        help_text = "Common values: Safety (0.45-0.60), Donut (0.70-0.80), Automatic (0.80-1.00)"
        ttk.Label(main_frame, text=help_text, font=('Arial', 8), foreground='gray').grid(row=3, column=0, columnspan=2,
                                                                                         pady=2)

        # Separator
        ttk.Separator(main_frame, orient='horizontal').grid(row=4, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=10)

        # Equipment correction section
        equipment_label = ttk.Label(main_frame, text="Equipment Corrections", font=('Arial', 12, 'bold'))
        equipment_label.grid(row=5, column=0, columnspan=2, pady=5)

        # Borehole diameter correction (CB)
        ttk.Label(main_frame, text="Borehole Diameter Correction (CB):").grid(row=6, column=0, sticky=tk.W, pady=5)
        self.borehole_corr_var = tk.StringVar(value="1.0")
        borehole_frame = ttk.Frame(main_frame)
        borehole_frame.grid(row=6, column=1, sticky=tk.W, pady=5)
        ttk.Entry(borehole_frame, textvariable=self.borehole_corr_var, width=10).pack(side=tk.LEFT)
        ttk.Button(borehole_frame, text="?", width=3, command=self.show_borehole_info).pack(side=tk.LEFT, padx=5)

        # Sampling method correction (CS)
        ttk.Label(main_frame, text="Sampler Correction (CS):").grid(row=7, column=0, sticky=tk.W, pady=5)
        self.sampler_corr_var = tk.StringVar(value="1.0")
        sampler_frame = ttk.Frame(main_frame)
        sampler_frame.grid(row=7, column=1, sticky=tk.W, pady=5)
        ttk.Entry(sampler_frame, textvariable=self.sampler_corr_var, width=10).pack(side=tk.LEFT)
        ttk.Button(sampler_frame, text="?", width=3, command=self.show_sampler_info).pack(side=tk.LEFT, padx=5)

        # Rod length correction (CR)
        ttk.Label(main_frame, text="Rod Length Correction (CR):").grid(row=8, column=0, sticky=tk.W, pady=5)
        self.rod_corr_var = tk.StringVar(value="1.0")
        rod_frame = ttk.Frame(main_frame)
        rod_frame.grid(row=8, column=1, sticky=tk.W, pady=5)
        ttk.Entry(rod_frame, textvariable=self.rod_corr_var, width=10).pack(side=tk.LEFT)
        ttk.Button(rod_frame, text="?", width=3, command=self.show_rod_info).pack(side=tk.LEFT, padx=5)

        # Equipment correction section
        equipment_label = ttk.Label(main_frame, text="Overburden Correction for N160", font=('Arial', 12, 'bold'))
        equipment_label.grid(row=9, column=0, columnspan=2, pady=5)

        # gwt input
        ttk.Label(main_frame, text="Groundwater table in feet:").grid(row=10, column=0, sticky=tk.W, pady=5)
        self.gwt_depth_var = tk.StringVar(value="2")
        ttk.Entry(main_frame, textvariable=self.gwt_depth_var, width=10).grid(row=10, column=1, sticky=(tk.W), pady=5)

        # Unit weight
        ttk.Label(main_frame, text="Unit Weight of Soil (pcf):").grid(row=11, column=0, sticky=tk.W, pady=5)
        self.unit_weight_var = tk.StringVar(value="120")
        ttk.Entry(main_frame, textvariable=self.unit_weight_var, width=10).grid(row=11, column=1, sticky=tk.W, pady=5)

        # Atmospheric pressure
        ttk.Label(main_frame, text="Atmospheric Pressure (psf):").grid(row=12, column=0, sticky=tk.W, pady=5)
        self.atm_pressure_var = tk.StringVar(value="2116")
        ttk.Entry(main_frame, textvariable=self.atm_pressure_var, width=10).grid(row=12, column=1, sticky=tk.W, pady=5)

        # Buttons frame
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=15, column=0, columnspan=2, pady=15)
        ttk.Button(button_frame, text="Calculate N60 & N160", command=self.calculate_n60).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Proceed to Modulus", command=self.open_modulus).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Back to Main Menu", command=self.back_to_main).pack(side=tk.LEFT, padx=5)

        self.results_text = scrolledtext.ScrolledText(main_frame, width=100, height=15, wrap=tk.WORD)
        self.results_text.grid(row=16, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)

        # Configure text tags for formatting
        self.results_text.tag_configure("header", font=('Courier', 9, 'bold'))
        self.results_text.tag_configure("data", font=('Courier', 9))

        main_frame.rowconfigure(16, weight=1)

    def back_to_main(self):
        self.root.destroy()
        main_menu()

    def show_borehole_info(self):
        """Display information about borehole diameter correction factors"""
        info = """Borehole Diameter Correction (CB):

Common values:
• 65-115 mm (2.5-4.5 in): CB = 1.00
• 150 mm (6 in): CB = 1.05
• 200 mm (8 in): CB = 1.15

Standard is typically 1.00 for normal borehole sizes."""
        messagebox.showinfo("Borehole Diameter Correction", info)

    def show_sampler_info(self):
        """Display information about sampler correction factors"""
        info = """Sampler Correction (CS):

Common values:
• Standard sampler: CS = 1.00
• Sampler without liner: CS = 1.10 to 1.30
• Sampler with liner: CS = 0.80 to 1.00

Standard split-spoon sampler with liner is typically 1.00."""
        messagebox.showinfo("Sampler Correction", info)

    def show_rod_info(self):
        """Display information about rod length correction factors"""
        info = """Rod Length Correction (CR):

Common values based on rod length:
• 3-4 m (10-13 ft): CR = 0.75
• 4-6 m (13-20 ft): CR = 0.85
• 6-10 m (20-33 ft): CR = 0.95
• 10-30 m (33-98 ft): CR = 1.00
• > 30 m (> 98 ft): CR = 1.00

Standard is 1.00 for rods longer than 10 m."""
        messagebox.showinfo("Rod Length Correction", info)

    def open_modulus(self):
        if not self.borings_data:
            messagebox.showerror("Error", "Please calculate N60 values first")
            return

        self.root.destroy()
        modulus_root = tk.Tk()
        DrainedModulus(modulus_root, self.borings_data)
        modulus_root.mainloop()

    def calculate_n60(self):
        if not self.borings_data:
            messagebox.showerror("Error", "Please add at least one boring")
            return

        try:
            hammer_eff = float(self.hammer_eff_var.get())
            unit_weight = float(self.unit_weight_var.get())
            atm_pressure = float(self.atm_pressure_var.get())
            borehole_corr = float(self.borehole_corr_var.get())
            sampler_corr = float(self.sampler_corr_var.get())
            rod_corr = float(self.rod_corr_var.get())
            gwt_depth = float(self.gwt_depth_var.get())

            if hammer_eff <= 0 or hammer_eff > 1:
                messagebox.showerror("Error", "Hammer efficiency must be between 0 and 1")
                return

            if unit_weight <= 0:
                messagebox.showerror("Error", "Unit weight must be greater than 0")
                return

            if atm_pressure <= 0:
                messagebox.showerror("Error", "Atmospheric pressure must be greater than 0")
                return

            self.results_text.delete(1.0, tk.END)

            for boring_id, data in self.borings_data.items():
                self.results_text.insert(tk.END, f"\n{'=' * 110}\n", "header")
                self.results_text.insert(tk.END, f"BORING: {boring_id}\n", "header")
                self.results_text.insert(tk.END, f"{'=' * 110}\n", "header")

                header = f"{'Depth':<8} {'N-fld':<8} {'σ\'v':<10} {'ER':<6} {'CB':<6} {'CS':<6} {'CR':<6} {'CN':<6} {'N60':<10} {'N160':<10}\n"
                self.results_text.insert(tk.END, header, "header")
                self.results_text.insert(tk.END, "-" * 110 + "\n", "header")

                for depth, n_field in zip(data['depths'], data['n_values']):
                    # Calculate effective overburden pressure
                    if gwt_depth >= depth:
                        eff_overburden = unit_weight * depth
                    else:
                        eff_overburden = (unit_weight * gwt_depth) + ((unit_weight-62.4) * (depth - gwt_depth))
                    # Energy ratio
                    ER = hammer_eff / 0.6

                    # Calculate N60 with all equipment corrections
                    n60 = convert_spt_to_n60(
                        n_field=n_field,
                        hammer_efficiency=hammer_eff,
                        borehole_diameter=borehole_corr,
                        sampler_correction=sampler_corr,
                        rod_length_correction=rod_corr
                    )

                    # Calculate N160 with overburden correction
                    n160 = convert_spt_to_n160(
                        n_field=n_field,
                        hammer_efficiency=hammer_eff,
                        overburden_pressure=eff_overburden,
                        atmospheric_pressure=atm_pressure,
                        borehole_diameter=borehole_corr,
                        sampler_correction=sampler_corr,
                        rod_length_correction=rod_corr
                    )

                    # Store N60 and N160 values back into borings_data
                    if 'n60_values' not in data:
                        data['n60_values'] = []
                    if 'n160_values' not in data:
                        data['n160_values'] = []

                    data['n60_values'].append(n60)
                    data['n160_values'].append(n160)

                    CN = (atm_pressure / eff_overburden) ** 0.5

                    line = f"{depth:<8.1f} {n_field:<8} {eff_overburden:<10.2f} {ER:<6.3f} {borehole_corr:<6.2f} {sampler_corr:<6.2f} {rod_corr:<6.2f} {CN:<6.3f} {n60:<10.2f} {n160:<10.2f}\n"
                    self.results_text.insert(tk.END, line, "data")
                    self.results_text.insert(tk.END, "\n")

            # Add summary
            self.results_text.insert(tk.END, f"\n{'=' * 110}\n", "header")
            self.results_text.insert(tk.END, "CORRECTION PARAMETERS USED:\n", "header")
            self.results_text.insert(tk.END, f"{'=' * 110}\n", "header")
            self.results_text.insert(tk.END, f"Energy Correction:\n", "data")
            self.results_text.insert(tk.END, f"  Hammer Efficiency: {hammer_eff:.2f} ({hammer_eff * 100:.0f}%)\n",
                                     "data")
            self.results_text.insert(tk.END, f"  Energy Ratio (ER): {hammer_eff / 0.6:.3f}\n", "data")
            self.results_text.insert(tk.END, f"\nEquipment Corrections:\n", "data")
            self.results_text.insert(tk.END, f"  Borehole Diameter (CB): {borehole_corr:.2f}\n", "data")
            self.results_text.insert(tk.END, f"  Sampler (CS): {sampler_corr:.2f}\n", "data")
            self.results_text.insert(tk.END, f"  Rod Length (CR): {rod_corr:.2f}\n", "data")
            self.results_text.insert(tk.END, f"\nOverburden Correction:\n", "data")
            self.results_text.insert(tk.END, f"  Unit Weight: {unit_weight:.1f} kN/m³\n", "data")
            self.results_text.insert(tk.END, f"  Atmospheric Pressure: {atm_pressure:.1f} kPa\n", "data")
            self.results_text.insert(tk.END, f"\nFormulas:\n", "data")
            self.results_text.insert(tk.END, f"  N60 = N-field × ER × CB × CS × CR\n", "data")
            self.results_text.insert(tk.END, f"  CN = (Pa / σ'v)^0.5\n", "data")
            self.results_text.insert(tk.END, f"  N160 = N60 × CN\n", "data")
            self.results_text.insert(tk.END, f"  σ'v = Unit Weight × Depth\n", "data")
            self.results_text.insert(tk.END, f"\nNote: Depth units in feet, σ'v in kPa\n", "data")

        except ValueError:
            messagebox.showerror("Error", "Invalid input values")

class SPT:

    def __init__(self, root):
        self.root = root
        self.root.title("Stratum Information")
        self.root.geometry("800x400")

        self.borings_data = {}

        # Create main frame
        main_frame = ttk.Frame(root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Configure grid weights
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)

        # Title
        title = ttk.Label(main_frame, text="Boring Information", font=('Arial', 16, 'bold'))
        title.grid(row=0, column=0, columnspan=2, pady=10)

        # Boring ID input
        ttk.Label(main_frame, text="Boring ID:").grid(row=1, column=0, sticky=tk.W, pady=5)
        self.boring_id_var = tk.StringVar(value="B-1")
        ttk.Entry(main_frame, textvariable=self.boring_id_var, width=20).grid(row=1, column=1, sticky=tk.W, pady=5)

        # N-values input
        ttk.Label(main_frame, text="N-values (comma-separated):").grid(row=2, column=0, sticky=tk.W, pady=5)
        self.n_values_var = tk.StringVar(value="10,29,40,17")
        ttk.Entry(main_frame, textvariable=self.n_values_var, width=40).grid(row=2, column=1, sticky=(tk.W, tk.E),
                                                                             pady=5)
        # Soil-type input
        ttk.Label(main_frame, text="Soil type (comma-separated, See note):").grid(row=3, column=0, sticky=tk.W, pady=5)
        self.soil_type_var = tk.StringVar(value="1,2,1,1")
        ttk.Entry(main_frame, textvariable=self.soil_type_var, width=40).grid(row=3, column=1, sticky=(tk.W, tk.E),
                                                                              pady=5)
        # Depths input
        ttk.Label(main_frame, text="Depths in feet (comma-separated):").grid(row=4, column=0, sticky=tk.W, pady=5)
        self.depths_var = tk.StringVar(value="1,3,5,7")
        ttk.Entry(main_frame, textvariable=self.depths_var, width=40).grid(row=4, column=1, sticky=(tk.W, tk.E), pady=5)

        # Add boring button
        ttk.Button(main_frame, text="Add Boring", command=self.add_boring).grid(row=6, column=0, columnspan=3, pady=10)

        # Soil type note
        ttk.Button(main_frame, text="Soil Type Note", command=self.soil_type_text_box).grid(row=3, column=2,
                                                                                            columnspan=1, pady=10)

        # Separator
        ttk.Separator(main_frame, orient='horizontal').grid(row=7, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=10)

        # Buttons frame
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=10, column=0, columnspan=2, pady=15)

        ttk.Button(button_frame, text="Clear All", command=self.clear_all).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Proceed to N60", command=self.open_n60).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Back to Main Menu", command=self.back_to_main).pack(side=tk.LEFT, padx=5)

        main_frame.rowconfigure(11, weight=1)

    def open_n60(self):
        if not self.borings_data:
            messagebox.showerror("Error", "Please add at least one boring before proceeding to N60 calculator")
            return

        self.root.destroy()
        n60_root = tk.Tk()
        N60(n60_root, self.borings_data)
        n60_root.mainloop()

    def back_to_main(self):
        self.root.destroy()
        main_menu()

    def add_boring(self):
        boring_id = self.boring_id_var.get().strip()
        n_values_str = self.n_values_var.get().strip()
        depths_str = self.depths_var.get().strip()
        soil_type_text = self.soil_type_var.get().strip()

        try:
            soil_type = [int(x.strip()) for x in soil_type_text.split(',')]
            for value in soil_type:
                if value < 1 or value > 6:
                    messagebox.showerror("Error", "Soil type must be between 1 and 6")
                    return
        except ValueError:
            messagebox.showerror("Error", "Invalid soil type value")
            return

        if not boring_id:
            messagebox.showerror("Error", "Please enter a Boring ID")
            return

        if not n_values_str or not depths_str:
            messagebox.showerror("Error", "Please enter both N-values and depths")
            return

        try:
            n_values = [int(x.strip()) for x in n_values_str.split(',')]
            depths = [float(x.strip()) for x in depths_str.split(',')]
            soil_types = [int(x.strip()) for x in soil_type_text.split(',')]

            if len(n_values) != len(depths) or len(n_values) != len(soil_types):
                messagebox.showerror("Error", "Number of values must match")
                return

            self.borings_data[boring_id] = {
                'n_values': n_values,
                'depths': depths,
                'soil types': soil_types
            }

            messagebox.showinfo("Success", f"Boring {boring_id} added with {len(n_values)} readings")

            # Clear input fields
            self.boring_id_var.set("")
            self.n_values_var.set("")
            self.depths_var.set("")
            self.soil_type_var.set("")

        except ValueError:
            messagebox.showerror("Error", "Invalid input. Please use numbers only, separated by commas")

    def clear_all(self):
        self.borings_data = {}
        self.boring_id_var.set("")
        self.n_values_var.set("")
        self.depths_var.set("")
        messagebox.showinfo("Cleared", "All data has been cleared")

    def soil_type_text_box(self):
        messagebox.showinfo("Soil Types",
                            "Type the number for the corresponding soil type:\n1: Sand (NC)\n2:Sand (saturated)\n3:Sand (OC)\n4:Gravelly Sand\n5:Clayey Sand\n6:Silts, sandy silts or clayey silt")

def main_menu():
    """Create the main menu window with navigation buttons"""
    root = tk.Tk()
    root.title("SPT Calculator - Main Menu")
    root.geometry("400x300")

    # Create main frame
    main_frame = ttk.Frame(root, padding="20")
    main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

    # Configure grid weights
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)
    main_frame.columnconfigure(0, weight=1)

    # Title
    title = ttk.Label(main_frame, text="SPT Calculator", font=('Arial', 20, 'bold'))
    title.grid(row=0, column=0, pady=30)

    # Subtitle
    subtitle = ttk.Label(main_frame, text="Select an option to continue:", font=('Arial', 12))
    subtitle.grid(row=1, column=0, pady=10)

    # Button frame for centering
    button_frame = ttk.Frame(main_frame)
    button_frame.grid(row=2, column=0, pady=20)

    # Subsurface Information button
    def open_subsurface():
        root.destroy()
        spt_root = tk.Tk()
        SPT(spt_root)
        spt_root.mainloop()

    subsurface_btn = ttk.Button(button_frame, text="Subsurface Information",
                                command=open_subsurface, width=25)
    subsurface_btn.pack(pady=10)

    # N60 Calculator button
    def open_n60():
        root.destroy()
        n60_root = tk.Tk()
        N60(n60_root)
        n60_root.mainloop()

    n60_btn = ttk.Button(button_frame, text="N60",
                         command=open_n60, width=25)
    n60_btn.pack(pady=10)

    # Drained Modulus Calculator button
    def open_modulus():
        root.destroy()
        modulus_root = tk.Tk()
        DrainedModulus(modulus_root)
        modulus_root.mainloop()

    modulus_btn = ttk.Button(button_frame, text="Drained Modulus Stiffness Calculator",
                             command=open_modulus, width=25)
    modulus_btn.pack(pady=10)

    # Exit button
    exit_btn = ttk.Button(button_frame, text="Exit",
                          command=root.destroy, width=25)
    exit_btn.pack(pady=10)

    root.mainloop()

class TerzaghiBearingCapacity:
    """
    A class to calculate ultimate bearing capacity using Terzaghi's bearing capacity equation.

    The general bearing capacity equation is:
    qu = c'*Nc*Fcs*Fcd*Fci + q*Nq*Fqs*Fqd*Fqi + 0.5*γ*B*Nγ*Fγs*Fγd*Fγi

    Where:
    - c' = cohesion
    - q = effective stress at foundation level (γ * D)
    - γ = unit weight of soil
    - B = width of foundation (diameter for circular foundation)
    - Nc, Nq, Nγ = bearing capacity factors
    - Fcs, Fqs, Fγs = shape factors
    - Fcd, Fqd, Fγd = depth factors
    - Fci, Fqi, Fγi = load inclination factors
    """

    def __init__(self, cohesion, friction_angle, unit_weight, foundation_depth,
                 foundation_width, gwt_depth, foundation_length=None, load_inclination=0):
        """
        Initialize the bearing capacity calculator.

        Parameters:
        -----------
        cohesion : float
            Cohesion of soil (c') in kPa or psf
        friction_angle : float
            Internal friction angle (φ') in degrees
        unit_weight : float
            Unit weight of soil (γ) in kN/m³ or pcf
        foundation_depth : float
            Depth of foundation (D) in m or ft
        foundation_width : float
            Width of foundation (B) in m or ft (diameter for circular)
        foundation_length : float, optional
            Length of foundation (L) in m or ft. If None, assumes square/circular foundation
        load_inclination : float, optional
            Inclination of load from vertical (β) in degrees. Default is 0 (vertical load)
        gwt_depth : flat
            Depth to the groundwater table
        """
        self.c = cohesion
        self.phi = friction_angle
        self.gamma = unit_weight
        self.D = foundation_depth
        self.B = foundation_width
        self.L = foundation_length if foundation_length else foundation_width
        self.beta = load_inclination
        self.gwt_depth = gwt_depth
        #TODO: Find a way to relate moist and saturated unit weight to relative density? Worth asking Myron
        self.gamma_eff = self.gamma_sat - 62.4
        self.gamma_sat = None

        # Convert angles to radians
        self.phi_rad = math.radians(friction_angle)
        self.beta_rad = math.radians(load_inclination)

    def groundwater_corrections(self):
        # Calculate effective stress at foundation level considering groundwater effects
        if self.gwt_depth <= self.D:
            self.q = (self.gamma * self.gwt_depth) + (self.gamma_eff * (self.D - self.gwt_depth))
            self.gamma = self.gamma_eff
        if self.gwt_depth <= (self.D + self.B):
            self.gamma = self.gamma_eff+((self.gwt_depth - self.D)/(self.B))*(self.gamma - self.gamma_eff)
    def calculate_bearing_capacity_factors(self):
        """
        Calculate Nc, Nq, and Nγ bearing capacity factors.

        Returns:
        --------
        tuple : (Nc, Nq, Nγ)
        """
        # Nq = tan²(45 + φ'/2) * e^(π * tan φ')
        Nq = (math.tan(math.radians(45) + self.phi_rad / 2) ** 2) * math.exp(math.pi * math.tan(self.phi_rad))

        # Nc = (Nq - 1) * cot φ'
        if self.phi > 0:
            Nc = (Nq - 1) * math.tan(self.phi_rad)
        else:
            Nc = 5.14  # For φ = 0 (purely cohesive soil)

        # Nγ = 2(Nq + 1) * tan φ'
        Ngamma = 2 * (Nq + 1) * math.tan(self.phi_rad)

        return Nc, Nq, Ngamma

    def calculate_shape_factors(self):
        """
        Calculate shape factors Fcs, Fqs, Fγs using DeBeer (1970) relationships.

        Returns:
        --------
        tuple : (Fcs, Fqs, Fγs)
        """

        Nc, Nq, _ = self.calculate_bearing_capacity_factors()

        # Fcs = 1 + (B/L) * (Nq/Nc)
        Fcs = 1 + (self.B / self.L) * (Nq / Nc)

        # Fqs = 1 + (B/L) * tan φ'
        Fqs = 1 + (self.B / self.L) * math.tan(self.phi_rad)

        # Fγs = 1 - 0.4 * (B/L)
        Fgammas = 1 - 0.4 * (self.B / self.L)

        return Fcs, Fqs, Fgammas
    def calculate_depth_factors(self):
        """
        Calculate depth factors Fcd, Fqd, Fγd using Hansen (1970) relationships.

        Returns:
        --------
        tuple : (Fcd, Fqd, Fγd)
        """
        Nc, Nq, _ = self.calculate_bearing_capacity_factors()
        D_B_ratio = self.D / self.B
        D_B_ratio_rad = math.radians(D_B_ratio)
        if D_B_ratio <= 1:
            if self.phi == 0:
                # For φ = 0
                Fcd = 1 + 0.4 * D_B_ratio
                Fqd = 1
                Fgammad = 1
            else:
                # For φ > 0
                Fqd=1+2*math.tan(self.phi_rad)*(1-math.sin(self.phi_rad)) ** 2 * D_B_ratio
                Fcd = Fqd - (1-Fqd)/(Nc*math.tan(self.phi_rad))
                Fgammad = 1
        else:
            if self.phi == 0:
                # For φ = 0
                Fcd = 1 + 0.4 * math.atan(D_B_ratio_rad)
                Fqd = 1
                Fgammad = 1
            else:
                # For φ > 0
                Fqd = 1 + 2 * math.tan(self.phi_rad) * (1 - math.sin(self.phi_rad)) ** 2 * math.atan(D_B_ratio_rad)
                Fcd = Fqd - (1 - Fqd) / (Nc * math.tan(self.phi_rad))
                Fgammad = 1

        return Fcd, Fqd, Fgammad

    def calculate_inclination_factors(self):
        """
        Calculate load inclination factors Fci, Fqi, Fγi using Meyerhof (1963) relationships.

        Returns:
        --------
        tuple : (Fci, Fqi, Fγi)
        """
        # Fci = Fqi = (1 - β°/90°)²
        Fci = (1 - self.beta / 90) ** 2
        Fqi = (1 - self.beta / 90) ** 2

        # Fγi = (1 - β/φ)²
        if self.phi > 0:
            Fgammai = (1 - self.beta_rad / self.phi_rad) ** 2
        else:
            Fgammai = 1  # For φ = 0

        return Fci, Fqi, Fgammai

    def calculate_ultimate_bearing_capacity(self):
        """
        Calculate the ultimate bearing capacity using the general equation.

        Returns:
        --------
        dict : Dictionary containing:
            - qu: ultimate bearing capacity
            - Nc, Nq, Nγ: bearing capacity factors
            - All shape, depth, and inclination factors
        """

        self.groundwater_corrections()
        
        #Calculate bearing capacity factors
        Nc, Nq, Ngamma = self.calculate_bearing_capacity_factors()

        # Calculate shape factors
        Fcs, Fqs, Fgammas = self.calculate_shape_factors()

        # Calculate depth factors
        Fcd, Fqd, Fgammad = self.calculate_depth_factors()

        # Calculate inclination factors
        Fci, Fqi, Fgammai = self.calculate_inclination_factors()

        # Calculate ultimate bearing capacity
        # qu = c'*Nc*Fcs*Fcd*Fci + q*Nq*Fqs*Fqd*Fqi + 0.5*γ*B*Nγ*Fγs*Fγd*Fγi
        term1 = self.c * Nc * Fcs * Fcd * Fci
        term2 = self.q * Nq * Fqs * Fqd * Fqi
        term3 = 0.5 * self.gamma * self.B * Ngamma * Fgammas * Fgammad * Fgammai

        qu = term1 + term2 + term3

        return {
            'qu': qu,
            'Nc': Nc,
            'Nq': Nq,
            'Ngamma': Ngamma,
            'Fcs': Fcs,
            'Fqs': Fqs,
            'Fgammas': Fgammas,
            'Fcd': Fcd,
            'Fqd': Fqd,
            'Fgammad': Fgammad,
            'Fci': Fci,
            'Fqi': Fqi,
            'Fgammai': Fgammai,
            'term1_cohesion': term1,
            'term2_surcharge': term2,
            'term3_self_weight': term3
        }

    def calculate_allowable_bearing_capacity(self, factor_of_safety=3.0):
        """
        Calculate the allowable bearing capacity.

        Parameters:
        -----------
        factor_of_safety : float, optional
            Factor of safety to apply. Default is 3.0

        Returns:
        --------
        float : Allowable bearing capacity
        """
        results = self.calculate_ultimate_bearing_capacity()
        qa = results['qu'] / factor_of_safety
        return qa

    def print_results(self, factor_of_safety=3.0):
        """
        Print a formatted summary of the bearing capacity calculation.

        Parameters:
        -----------
        factor_of_safety : float, optional
            Factor of safety to apply. Default is 3.0
        """
        results = self.calculate_ultimate_bearing_capacity()
        qa = self.calculate_allowable_bearing_capacity(factor_of_safety)

        print("=" * 70)
        print("TERZAGHI BEARING CAPACITY CALCULATION")
        print("=" * 70)
        print(f"\nInput Parameters:")
        print(f"  Cohesion (c'):              {self.c:.2f}")
        print(f"  Friction Angle (φ'):        {self.phi:.2f}°")
        print(f"  Unit Weight (γ):            {self.gamma:.2f}")
        print(f"  Foundation Depth (D):       {self.D:.2f}")
        print(f"  Foundation Width (B):       {self.B:.2f}")
        print(f"  Foundation Length (L):      {self.L:.2f}")
        print(f"  Load Inclination (β):       {self.beta:.2f}°")
        print(f"  Effective Stress (q):       {self.q:.2f}")

        print(f"\nBearing Capacity Factors:")
        print(f"  Nc = {results['Nc']:.3f}")
        print(f"  Nq = {results['Nq']:.3f}")
        print(f"  Nγ = {results['Ngamma']:.3f}")

        print(f"\nShape Factors:")
        print(f"  Fcs = {results['Fcs']:.3f}")
        print(f"  Fqs = {results['Fqs']:.3f}")
        print(f"  Fγs = {results['Fgammas']:.3f}")

        print(f"\nDepth Factors:")
        print(f"  Fcd = {results['Fcd']:.3f}")
        print(f"  Fqd = {results['Fqd']:.3f}")
        print(f"  Fγd = {results['Fgammad']:.3f}")

        print(f"\nInclination Factors:")
        print(f"  Fci = {results['Fci']:.3f}")
        print(f"  Fqi = {results['Fqi']:.3f}")
        print(f"  Fγi = {results['Fgammai']:.3f}")

        print(f"\nBearing Capacity Components:")
        print(f"  Cohesion term:      {results['term1_cohesion']:.2f}")
        print(f"  Surcharge term:     {results['term2_surcharge']:.2f}")
        print(f"  Self-weight term:   {results['term3_self_weight']:.2f}")

        print(f"\nResults:")
        print(f"  Ultimate Bearing Capacity (qu): {results['qu']:.2f}")
        print(f"  Allowable Bearing Capacity (qa) with FS={factor_of_safety}: {qa:.2f}")
        print("=" * 70)

class Stratums_and_SoilProps:
    # Typical soil property correlations (Imperial units)
    # Based on numerical soil type codes
    # TODO: check these correlations
    SOIL_PROPERTIES = {
        1: {  # Sand (NC - Normally Consolidated)
            'name': 'Sand (NC)',
            'unit_weight': lambda n: 100 + 1.0 * n,  # Loose to dense (pcf)
            'phi': lambda n: 28 + 0.4 * n,  # Common correlation for NC sand
            'cohesion': lambda n: 0,  # Cohesionless
            'delta': lambda n: (28 + 0.4 * n) / 2
        },
        2: {  # Sand (saturated)
            'name': 'Sand (saturated)',
            'unit_weight': lambda n: 115 + 0.5 * n,  # Saturated unit weight (pcf)
            'phi': lambda n: 28 + 0.35 * n,  # Slightly reduced for saturated conditions
            'cohesion': lambda n: 0,  # Cohesionless
            'delta': lambda n: (28 + 0.35 * n) / 2
        },
        3: {  # Sand (OC - Overconsolidated)
            'name': 'Sand (OC)',
            'unit_weight': lambda n: 105 + 1.2 * n,  # Denser for OC (pcf)
            'phi': lambda n: 30 + 0.45 * n,  # Higher friction angle for OC
            'cohesion': lambda n: 0,  # Cohesionless
            'delta': lambda n: (30 + 0.45 * n) / 2
        },
        4: {  # Gravelly Sand
            'name': 'Gravelly Sand',
            'unit_weight': lambda n: 115 + 1.0 * n,  # Higher density (pcf)
            'phi': lambda n: 32 + 0.4 * n,  # Higher friction angle
            'cohesion': lambda n: 0,  # Cohesionless
            'delta': lambda n: (32 + 0.4 * n) / 2
        },
        5: {  # Clayey Sand
            'name': 'Clayey Sand',
            'unit_weight': lambda n: 105 + 0.8 * n,  # pcf
            'phi': lambda n: 26 + 0.35 * n,  # Reduced friction angle
            'cohesion': lambda n: 0.02 + 0.003 * n,  # Small cohesion (ksf)
            'delta': lambda n: (26 + 0.35 * n) / 2
        },
        6: {  # Silts, sandy silts or clayey silt
            'name': 'Silts, sandy silts or clayey silt',
            'unit_weight': lambda n: 100 + 0.7 * n,  # pcf
            'phi': lambda n: 24 + 0.3 * n,  # Lower friction angle
            'cohesion': lambda n: 0.025 + 0.004 * n,  # Moderate cohesion (ksf)
            'delta': lambda n: (24 + 0.3 * n) / 2
        }
    }

    def __init__(self, thickness: float = None, unit_weight: float = None,
                 phi: float = None, cohesion: float = 0, delta: Optional[float] = None,
                 n_value: float = None, soil_type: int = None):
        """
        Initialize a soil layer either with direct properties or from N-value correlations.

        Args:
            thickness: Layer thickness (ft)
            unit_weight: Unit weight (pcf - pounds per cubic foot)
            phi: Friction angle (degrees)
            cohesion: Cohesion (ksf - kips per square foot)
            delta: Wall friction angle (degrees), defaults to phi/2
            n_value: SPT N-value for property correlation
            soil_type: Soil type code (1-6):
                1: Sand (NC)
                2: Sand (saturated)
                3: Sand (OC)
                4: Gravelly Sand
                5: Clayey Sand
                6: Silts, sandy silts or clayey silt
        """
        self.thickness = thickness

        # If N-value and soil type provided, use correlations
        if n_value is not None and soil_type is not None:
            # Validate soil type code
            if soil_type not in self.SOIL_PROPERTIES:
                raise ValueError(f"Invalid soil type code: {soil_type}. Must be 1-6.")

            props = self.SOIL_PROPERTIES[soil_type]

            self.soil_type_name = props['name']
            self.unit_weight = unit_weight if unit_weight is not None else props['unit_weight'](n_value)
            self.phi = phi if phi is not None else props['phi'](n_value)
            self.cohesion = cohesion if cohesion is not None else props['cohesion'](n_value)
            self.delta = delta if delta is not None else props['delta'](n_value)
        else:
            # Use provided values or defaults
            self.soil_type_name = 'Unspecified'
            self.unit_weight = unit_weight if unit_weight is not None else 115.0  # pcf
            self.phi = phi if phi is not None else 30.0
            self.cohesion = cohesion
            self.delta = delta if delta is not None else self.phi / 2

    @classmethod
    def from_boring_data(cls, boring_data: dict, boring_id: str) -> List['Stratums_and_SoilProps']:
        """
        Create soil layers from boring data dictionary.

        Args:
            boring_data: Dictionary with structure:
                {boring_id: {
                    'n_values': [list of N values],
                    'depths': [list of depths in feet],
                    'soil_types': [list of soil type codes 1-6]
                }}
            boring_id: ID of the boring to use

        Returns:
            List of SoilLayer objects
        """
        if boring_id not in boring_data:
            raise ValueError(f"Boring ID '{boring_id}' not found in boring data")

        data = boring_data[boring_id]
        n_values = data['n_values']
        depths = data['depths']
        soil_types = data['soil_types']

        # Validate data
        if not (len(n_values) == len(depths) == len(soil_types)):
            raise ValueError("n_values, depths, and soil_types must have same length")

        layers = []

        for i in range(len(depths)):
            # Calculate thickness
            if i == 0:
                thickness = depths[i]
            else:
                thickness = depths[i] - depths[i - 1]

            # Create layer using N-value correlation
            layer = cls(
                thickness=thickness,
                n_value=n_values[i],
                soil_type=soil_types[i]
            )

            layers.append(layer)

        return layers

    def __repr__(self):
        return (f"SoilLayer(type={self.soil_type_name}, thickness={self.thickness:.2f} ft, "
                f"γ={self.unit_weight:.1f} pcf, φ={self.phi:.1f}°, "
                f"c={self.cohesion:.3f} ksf)")


class LineLoad:
    def __init__(self, magnitude: float, distance: float):
        """
        Line load (point load).

        Args:
            magnitude: Load magnitude (kN/m)
            distance: Horizontal distance from wall (m)
        """
        self.magnitude = magnitude
        self.distance = distance

class DistributedLoad:
    def __init__(self, magnitude: float, start_distance: float, end_distance: float):
        """
        Distributed surcharge load.

        Args:
            magnitude: Load magnitude (kPa)
            start_distance: Start distance from wall (m)
            end_distance: End distance from wall (m)
        """
        self.magnitude = magnitude
        self.start_distance = start_distance
        self.end_distance = end_distance

class Strut:
    def __init__(self, depth: float, spacing: float = 1.0):
        """
        Strut/brace support.

        Args:
            depth: Depth from top of wall (m)
            spacing: Tributary spacing (m)
        """
        self.depth = depth
        self.spacing = spacing

class LateralEarthPressure:
    def __init__(self, soil_layers: List[Stratums_and_SoilProps], wall_height: float,
                 theory: str = 'rankine', line_loads: List[LineLoad] = None,
                 distributed_loads: List[DistributedLoad] = None,
                 struts: List[Strut] = None, wall_type: str = 'retaining'):
        """
        Initialize lateral earth pressure calculator.

        Args:
            soil_layers: List of soil layers from top to bottom
            wall_height: Total wall height (m)
            theory: 'rankine', 'coulomb', 'fhwa', or 'peck'
            line_loads: List of line loads
            distributed_loads: List of distributed loads
            struts: List of struts/braces
            wall_type: 'retaining' or 'braced'
        """
        self.soil_layers = soil_layers
        self.wall_height = wall_height
        self.theory = theory.lower()
        self.line_loads = line_loads if line_loads else []
        self.distributed_loads = distributed_loads if distributed_loads else []
        self.struts = struts if struts else []
        self.wall_type = wall_type
        self.excavation_depth = wall_height

    def calculate_rankine_ka(self, phi: float, beta: float = 0) -> float:
        """Calculate Rankine active earth pressure coefficient."""
        phi_rad = np.radians(phi)
        beta_rad = np.radians(beta)

        cos_beta = np.cos(beta_rad)
        cos_phi = np.cos(phi_rad)

        numerator = cos_beta * (cos_beta - np.sqrt(cos_beta ** 2 - cos_phi ** 2))
        denominator = cos_beta + np.sqrt(cos_beta ** 2 - cos_phi ** 2)

        return numerator / denominator

    def calculate_rankine_kp(self, phi: float, beta: float = 0) -> float:
        """Calculate Rankine passive earth pressure coefficient."""
        phi_rad = np.radians(phi)
        beta_rad = np.radians(beta)

        cos_beta = np.cos(beta_rad)
        cos_phi = np.cos(phi_rad)

        numerator = cos_beta * (cos_beta + np.sqrt(cos_beta ** 2 - cos_phi ** 2))
        denominator = cos_beta - np.sqrt(cos_beta ** 2 - cos_phi ** 2)

        return numerator / denominator

    #TODO: The form of the coulomb formulas is correct but the terms are in the wrong spot. Fix and rename variables to convention
    def calculate_coulomb_ka(self, phi: float, delta: float,
                             beta: float = 0, alpha: float = 90) -> float:
        """Calculate Coulomb active earth pressure coefficient."""
        phi_rad = np.radians(phi)
        delta_rad = np.radians(delta)
        beta_rad = np.radians(beta)
        alpha_rad = np.radians(alpha)

        sin_phi_plus_delta = np.sin(phi_rad + delta_rad)
        sin_phi_minus_beta = np.sin(phi_rad - beta_rad)
        sin_alpha_plus_delta = np.sin(alpha_rad + delta_rad)
        sin_alpha_minus_beta = np.sin(alpha_rad - beta_rad)

        sqrt_term = np.sqrt(sin_phi_plus_delta * sin_phi_minus_beta /
                            (sin_alpha_plus_delta * sin_alpha_minus_beta))

        numerator = np.sin(alpha_rad + phi_rad) ** 2
        denominator = (np.sin(alpha_rad) ** 2 * sin_alpha_minus_beta *
                       (1 + sqrt_term) ** 2)

        return numerator / denominator

    def calculate_coulomb_kp(self, phi: float, delta: float,
                             beta: float = 0, alpha: float = 90) -> float:
        """Calculate Coulomb passive earth pressure coefficient."""
        phi_rad = np.radians(phi)
        delta_rad = np.radians(delta)
        beta_rad = np.radians(beta)
        alpha_rad = np.radians(alpha)

        sin_phi_minus_delta = np.sin(phi_rad - delta_rad)
        sin_phi_plus_beta = np.sin(phi_rad + beta_rad)
        sin_alpha_minus_delta = np.sin(alpha_rad - delta_rad)
        sin_alpha_plus_beta = np.sin(alpha_rad + beta_rad)

        sqrt_term = np.sqrt(sin_phi_minus_delta * sin_phi_plus_beta /
                            (sin_alpha_minus_delta * sin_alpha_plus_beta))

        numerator = np.sin(alpha_rad - phi_rad) ** 2
        denominator = (np.sin(alpha_rad) ** 2 * sin_alpha_plus_beta *
                       (1 - sqrt_term) ** 2)

        return numerator / denominator

    def calculate_fhwa_envelope(self, depth: float, soil_type: str = 'sand') -> float:
        """Calculate FHWA apparent earth pressure envelope."""
        H = self.excavation_depth

        if soil_type == 'sand':
            gamma = self.get_average_unit_weight()
            phi = self.get_average_friction_angle()
            Ka = self.calculate_rankine_ka(phi)

            if depth <= 0.25 * H:
                pressure = 0.65 * Ka * gamma * H * (depth / (0.25 * H))
            else:
                pressure = 0.65 * Ka * gamma * H
        else:  # clay
            gamma = self.get_average_unit_weight()
            c = self.get_average_cohesion()
            pressure = max(0, gamma * H - 4 * c)

        return pressure

    def calculate_peck_envelope(self, depth: float, soil_type: str = 'sand') -> float:
        """Calculate Peck's apparent earth pressure envelope."""
        H = self.excavation_depth
        gamma = self.get_average_unit_weight()

        if soil_type == 'sand':
            Ka = 0.65 * self.calculate_rankine_ka(self.get_average_friction_angle())
            pressure = Ka * gamma * H
        else:  # clay
            c = self.get_average_cohesion()
            m = min(1.0, 4 / (self.excavation_depth / (2 * c / gamma)))
            pressure = gamma * H * m

        return pressure

    def get_average_unit_weight(self) -> float:
        """Get weighted average unit weight."""
        if not self.soil_layers:
            return 18.0
        total_depth = sum(layer.thickness for layer in self.soil_layers)
        weighted_sum = sum(layer.unit_weight * layer.thickness
                           for layer in self.soil_layers)
        return weighted_sum / total_depth

    def get_average_friction_angle(self) -> float:
        """Get weighted average friction angle."""
        if not self.soil_layers:
            return 30.0
        total_depth = sum(layer.thickness for layer in self.soil_layers)
        weighted_sum = sum(layer.phi * layer.thickness
                           for layer in self.soil_layers)
        return weighted_sum / total_depth

    def get_average_cohesion(self) -> float:
        """Get weighted average cohesion."""
        if not self.soil_layers:
            return 0.0
        total_depth = sum(layer.thickness for layer in self.soil_layers)
        weighted_sum = sum(layer.cohesion * layer.thickness
                           for layer in self.soil_layers)
        return weighted_sum / total_depth

    def boussinesq_line_load(self, load: float, distance: float, depth: float) -> float:
        """Calculate horizontal stress from line load using Boussinesq theory."""
        if depth <= 0:
            return 0

        Q = load
        x = distance
        z = depth

        r_squared = x ** 2 + z ** 2
        sigma_h = (2 * Q * x ** 2 * z) / (np.pi * r_squared ** 2)

        return sigma_h

    def boussinesq_distributed_load(self, load: float, start_dist: float,
                                    end_dist: float, depth: float) -> float:
        """Calculate horizontal stress from distributed load using Boussinesq theory."""
        if depth <= 0:
            return 0

        q = load
        z = depth

        beta1 = np.arctan2(start_dist, z)
        beta2 = np.arctan2(end_dist, z)

        sigma_h = (q / np.pi) * (beta2 - beta1 +
                                 0.5 * (np.sin(2 * beta2) - np.sin(2 * beta1)))

        return sigma_h

    def calculate_pressure_profile(self, side: str = 'active',
                                   num_points: int = 100) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate lateral earth pressure profile.

        Args:
            side: 'active' or 'passive'
            num_points: Number of points in profile

        Returns:
            Tuple of (depths, pressures) arrays
        """
        depths = np.linspace(0, self.wall_height, num_points)
        pressures = np.zeros(num_points)

        for i, depth in enumerate(depths):
            pressure = 0

            if self.theory == 'fhwa':
                pressure = self.calculate_fhwa_envelope(depth, 'sand')
            elif self.theory == 'peck':
                pressure = self.calculate_peck_envelope(depth, 'sand')
            else:
                # Rankine or Coulomb
                cumulative_depth = 0
                sigma_v = 0

                for layer in self.soil_layers:
                    if depth > cumulative_depth and depth <= cumulative_depth + layer.thickness:
                        depth_in_layer = depth - cumulative_depth
                        sigma_v += layer.unit_weight * depth_in_layer

                        if side == 'active':
                            if self.theory == 'coulomb':
                                K = self.calculate_coulomb_ka(layer.phi, layer.delta)
                            else:
                                K = self.calculate_rankine_ka(layer.phi)
                        else:  # passive
                            if self.theory == 'coulomb':
                                K = self.calculate_coulomb_kp(layer.phi, layer.delta)
                            else:
                                K = self.calculate_rankine_kp(layer.phi)

                        pressure = K * sigma_v - 2 * layer.cohesion * np.sqrt(K)
                        break
                    elif depth > cumulative_depth + layer.thickness:
                        sigma_v += layer.unit_weight * layer.thickness
                        cumulative_depth += layer.thickness

            # Add surcharge from line loads
            for line_load in self.line_loads:
                pressure += self.boussinesq_line_load(
                    line_load.magnitude, line_load.distance, depth)

            # Add surcharge from distributed loads
            for dist_load in self.distributed_loads:
                pressure += self.boussinesq_distributed_load(
                    dist_load.magnitude, dist_load.start_distance,
                    dist_load.end_distance, depth)

            pressures[i] = max(0, pressure)

        return depths, pressures

    def calculate_strut_reactions(self, depths: np.ndarray,
                                  pressures: np.ndarray) -> List[Dict]:
        """
        Calculate reactions at strut locations.

        Args:
            depths: Depth array from pressure profile
            pressures: Pressure array from pressure profile

        Returns:
            List of dicts with 'depth' and 'reaction' keys
        """
        reactions = []

        sorted_struts = sorted(self.struts, key=lambda s: s.depth)

        for i, strut in enumerate(sorted_struts):
            strut_depth = strut.depth
            spacing = strut.spacing

            # Interpolate pressure at strut depth
            pressure_at_strut = np.interp(strut_depth, depths, pressures)

            # Calculate tributary height
            if i == 0:
                if len(sorted_struts) > 1:
                    tributary_height = ((sorted_struts[1].depth - strut_depth) / 2 +
                                        strut_depth)
                else:
                    tributary_height = self.wall_height
            elif i == len(sorted_struts) - 1:
                tributary_height = ((strut_depth - sorted_struts[i - 1].depth) / 2 +
                                    (self.wall_height - strut_depth))
            else:
                tributary_height = ((strut_depth - sorted_struts[i - 1].depth) / 2 +
                                    (sorted_struts[i + 1].depth - strut_depth) / 2)

            # Calculate average pressure over tributary area
            start_depth = max(0, strut_depth - tributary_height / 2)
            end_depth = min(self.wall_height, strut_depth + tributary_height / 2)

            mask = (depths >= start_depth) & (depths <= end_depth)
            avg_pressure = np.mean(pressures[mask]) if np.any(mask) else 0

            reaction = avg_pressure * tributary_height * spacing

            reactions.append({
                'depth': strut_depth,
                'reaction': reaction
            })

        return reactions

    def calculate_shear_moment(self, depths: np.ndarray, pressures: np.ndarray,
                               strut_reactions: List[Dict]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate shear and moment diagrams.

        Args:
            depths: Depth array from pressure profile
            pressures: Pressure array from pressure profile
            strut_reactions: List of strut reactions

        Returns:
            Tuple of (shear, moment) arrays
        """
        num_points = len(depths)
        shear = np.zeros(num_points)
        moment = np.zeros(num_points)

        dz = self.wall_height / (num_points - 1)

        for i in range(num_points):
            depth = depths[i]
            pressure = pressures[i]

            # Check for strut at this depth
            strut_at_depth = None
            for strut_reaction in strut_reactions:
                if abs(strut_reaction['depth'] - depth) < dz / 2:
                    strut_at_depth = strut_reaction
                    break

            # Calculate incremental load
            if i > 0:
                avg_pressure = (pressures[i - 1] + pressure) / 2
                incremental_load = avg_pressure * dz

                shear[i] = shear[i - 1] + incremental_load
                moment[i] = moment[i - 1] + shear[i - 1] * dz + incremental_load * dz / 2

            # Apply strut reaction
            if strut_at_depth:
                shear[i] -= strut_at_depth['reaction']

        return shear, moment

if __name__ == "__main__":
    main_menu()