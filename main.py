import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext
import math

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

def convert_spt_to_n160(n_field, hammer_efficiency=0.6, overburden_pressure=100, atmospheric_pressure=100.0,
                        borehole_diameter=1.0, sampler_correction=1.0, rod_length_correction=1.0):
    """Convert SPT field N-value to N160 (N1)60 using energy, equipment, and overburden correction."""
    # First calculate N60 with all equipment corrections
    N60 = convert_spt_to_n60(n_field, hammer_efficiency, borehole_diameter, sampler_correction, rod_length_correction)

    # Overburden correction factor CN
    # CN = (Pa / σ'v)^0.5 where Pa is atmospheric pressure and σ'v is effective overburden
    CN = (atmospheric_pressure / overburden_pressure) ** 0.5

    # N160 = N60 * CN
    N160 = N60 * CN
    return N160

#todo: Review this matfoundation calculator
class MatFoundationSettlement:
    """
    Calculate settlement under a mat foundation using Boussinesq's theory
    and geotechnical soil parameters.
    """

    def __init__(self, width, length, depth, borings_data, E_data):
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
        self.E_data = E_data

        # Storage for loading conditions
        self.point_load = 0  # lbs
        self.uniform_load = 0  # psf

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

    def _boussinesq_uniform_load(self, q, B, L, z):
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
        m = L / z
        n = B / z
        #todo: go into bowles 1996 and confirm the below equation
        #Influence factor I for rectangular area (Fadum's chart approximation)
        term1 = 2 * m * n * math.sqrt(m ** 2 + n ** 2 + 1)
        term2 = (m ** 2 + n ** 2 + 2) * (m ** 2 + n ** 2 + 1)
        term3 = m ** 2 + n ** 2

        I = (1 / (4 * math.pi)) * (term1 / term2 + math.atan(term1 / (term3 * math.sqrt(m ** 2 + n ** 2 + 1))))

        sigma_z_pl = self.point_load / (self.length * self.width)
        sigma_z_tot = (q+sigma_z_pl) * I

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

        depths = self.borings_data[boring_id]['depth']
        stress_profile = []

        for depth in depths:
            # Adjust depth to account for foundation depth
            z = depth - self.depth

            if z <= 0:
                stress_profile.append(0)
                continue

            # Calculate stress from uniform load
            sigma_uniform = self._boussinesq_uniform_load(
                self.uniform_load, self.width, self.length, z
            )

            # Calculate stress from point load
            sigma_point = self.point_load / (self.length * self.width)

            # Total stress
            total_stress = sigma_uniform + sigma_point
            stress_profile.append(total_stress)

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
        if boring_id not in self.E_data:
            raise ValueError(f"Boring ID '{boring_id}' not found in E_data")

        depths = self.borings_data[boring_id]['depth']
        n60_values = self.borings_data[boring_id]['n60']
        E_values = self.E_data[boring_id]['E']
        stress_profile = self.calculate_stress_under_footing(boring_id)

        settlement_profile = []

        for i, (depth, stress, n60, E) in enumerate(zip(depths, stress_profile, n60_values, E_values)):

            # Determine layer thickness (distance to next depth or assumed increment)
            if i < len(depths) - 1:
                layer_thickness = depths[i + 1] - depth
            else:
                layer_thickness = 5  # Assume 5 ft for last layer

            # PLACEHOLDER: Replace with actual settlement formula
            settlement = (stress * layer_thickness) / E

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
        Calculate E (psf) based on soil type and N60 value.
        Each soil type uses a different formula.
        """
        if soil_type == 1:  # Sand (NC)
            E = 15 + 10*n60
            return E

        elif soil_type == 2:  # Sand (saturated)
            E = 15 * 5*n60
            return E

        elif soil_type == 3:  # Sand (OC)
            E = (40000 + 1050*n60)/50
            return E

        elif soil_type == 4:  # Gravelly Sand
            if n60 > 15:
                E = ((n60 + 6) * 600 + 2000) / 50
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
        #todo: save the values of E into something that can be used later on for calcs

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
        self.gwt_depth_var = tk.StringVar(value="100")
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

            n60_list = []
            n160_list = []

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
        self.boring_id_var = tk.StringVar()
        ttk.Entry(main_frame, textvariable=self.boring_id_var, width=20).grid(row=1, column=1, sticky=tk.W, pady=5)

        # N-values input
        ttk.Label(main_frame, text="N-values (comma-separated):").grid(row=2, column=0, sticky=tk.W, pady=5)
        self.n_values_var = tk.StringVar()
        ttk.Entry(main_frame, textvariable=self.n_values_var, width=40).grid(row=2, column=1, sticky=(tk.W, tk.E),
                                                                             pady=5)
        # Soil-type input
        ttk.Label(main_frame, text="Soil type (comma-separated, See note):").grid(row=3, column=0, sticky=tk.W, pady=5)
        self.soil_type_var = tk.StringVar()
        ttk.Entry(main_frame, textvariable=self.soil_type_var, width=40).grid(row=3, column=1, sticky=(tk.W, tk.E),
                                                                              pady=5)
        # Depths input
        ttk.Label(main_frame, text="Depths in feet (comma-separated):").grid(row=4, column=0, sticky=tk.W, pady=5)
        self.depths_var = tk.StringVar()
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

if __name__ == "__main__":
    main_menu()