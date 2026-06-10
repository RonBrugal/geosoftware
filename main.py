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
    def __init__(self, borings_data=None):
        self.borings_data = borings_data
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

    def run(self):
        print("\n" + "=" * 70)
        print("Drained Modulus Stiffness Calculator")
        print("=" * 70)
        print("Calculates E (kPa) using soil-specific formulas with N60 as input")

        while True:
            print("\nOptions:")
            print("  1. Calculate E")
            print("  2. View Formulas")
            print("  3. Proceed to Footing Settlement on Sand")
            print("  4. Back to Main Menu")
            choice = input("Enter choice: ").strip()

            if choice == '1':
                self.calculate_modulus()
            elif choice == '2':
                self.show_formulas()
            elif choice == '3':
                FootingSettlementSand(self.borings_data).run()
                return
            elif choice == '4':
                return
            else:
                print("Invalid choice")

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
        print("\nFormulas for Calculating E (kPa) from N60:")
        print("=" * 70)
        for soil_type, data in sorted(self.soil_formulas.items()):
            print(f"\nSoil Type {soil_type}: {data['name']}")
            print(f"Formula: {data['formula']}")
            print(f"Description: {data['description']}")
            print("-" * 70)

    def calculate_modulus(self):
        if not self.borings_data:
            print("Error: Please add boring data with soil types first")
            return

        for boring_id, data in self.borings_data.items():
            if 'n60_values' not in data:
                print(f"Error: No N60 values found for {boring_id}. Run N60 calculation first.")
                continue
            print(f"\n{'=' * 85}")
            print(f"BORING: {boring_id}")
            print(f"{'=' * 85}")
            print(f"{'Depth (ft)':<12} {'Soil Type':<12} {'Soil Name':<25} {'N60':<12} {'E (kPa)':<12}")
            print("-" * 85)

            for depth, soil_type, n60 in zip(data['depths'], data['soil types'], data['n60_values']):
                if soil_type in self.soil_formulas:
                    soil_name = self.soil_formulas[soil_type]['name']
                    E = self.calculate_e_from_formula(soil_type, n60) if n60 is not None else None

                    if 'E' not in data:
                        data['E'] = []
                    data['E'].append(E)

                    if E is not None:
                        print(f"{depth:<12.1f} {soil_type:<12} {soil_name:<25} {n60:<12.2f} {E:<12.2f}")
                    else:
                        print(f"{depth:<12.1f} {soil_type or 'N/A':<12} {soil_name:<25} {'N/A':<12} {'N/A':<12}")
                else:
                    print(f"{depth:<12.1f} {soil_type:<12} {'Unknown':<25} {n60:<12.2f} {'N/A':<12}")
            print()

        print(f"\n{'=' * 85}")
        print("CALCULATION NOTES:")
        print(f"{'=' * 85}")
        print("Each soil type uses a different formula to calculate E from N60")
        print("Select 'View Formulas' to see the formulas for each soil type")
        mat_calcs = MatFoundationSettlement(width=5, length=5, depth=2, borings_data=self.borings_data)
        for boring_id in self.borings_data.keys():
            print(mat_calcs.calculate_settlement(boring_id=boring_id))

class N60:
    def __init__(self, borings_data=None):
        self.borings_data = borings_data

    def run(self):
        print("\n" + "=" * 70)
        print("N60 and N160 Calculator")
        print("=" * 70)

        while True:
            print("\nOptions:")
            print("  1. Calculate N60 & N160")
            print("  2. View Borehole Diameter Correction (CB) info")
            print("  3. View Sampler Correction (CS) info")
            print("  4. View Rod Length Correction (CR) info")
            print("  5. Proceed to Drained Modulus")
            print("  6. Proceed to Liquefaction Analysis")
            print("  7. Back to Main Menu")
            choice = input("Enter choice: ").strip()

            if choice == '1':
                self.calculate_n60()
            elif choice == '2':
                self.show_borehole_info()
            elif choice == '3':
                self.show_sampler_info()
            elif choice == '4':
                self.show_rod_info()
            elif choice == '5':
                self.open_modulus()
                return
            elif choice == '6':
                self.open_liquefaction()
                return
            elif choice == '7':
                return
            else:
                print("Invalid choice")

    def show_borehole_info(self):
        print("""
Borehole Diameter Correction (CB):

Common values:
  - 65-115 mm (2.5-4.5 in): CB = 1.00
  - 150 mm (6 in): CB = 1.05
  - 200 mm (8 in): CB = 1.15

Standard is typically 1.00 for normal borehole sizes.
""")

    def show_sampler_info(self):
        print("""
Sampler Correction (CS):

Common values:
  - Standard sampler: CS = 1.00
  - Sampler without liner: CS = 1.10 to 1.30
  - Sampler with liner: CS = 0.80 to 1.00

Standard split-spoon sampler with liner is typically 1.00.
""")

    def show_rod_info(self):
        print("""
Rod Length Correction (CR):

Common values based on rod length:
  - 3-4 m (10-13 ft): CR = 0.75
  - 4-6 m (13-20 ft): CR = 0.85
  - 6-10 m (20-33 ft): CR = 0.95
  - 10-30 m (33-98 ft): CR = 1.00
  - > 30 m (> 98 ft): CR = 1.00

Standard is 1.00 for rods longer than 10 m.
""")

    def open_modulus(self):
        if not self.borings_data:
            print("Error: Please calculate N60 values first")
            return
        DrainedModulus(self.borings_data).run()

    def open_liquefaction(self):
        if not self.borings_data:
            print("Error: Please calculate N60 values first")
            return
        Liquefaction(self.borings_data, magnitude=7.5,
                     peak_acceleration_g=0.20, gwt_depth=5.0).run()

    def _prompt_float(self, label, default):
        value = input(f"{label} [{default}]: ").strip()
        return float(value) if value else float(default)

    def calculate_n60(self):
        if not self.borings_data:
            print("Error: Please add at least one boring")
            return

        try:
            print("\nEnter correction parameters (press Enter to accept default):")
            print("Hammer Efficiency common values: Safety (0.45-0.60), Donut (0.70-0.80), Automatic (0.80-1.00)")
            hammer_eff = self._prompt_float("Hammer Efficiency (decimal)", "0.80")
            borehole_corr = self._prompt_float("Borehole Diameter Correction (CB)", "1.0")
            sampler_corr = self._prompt_float("Sampler Correction (CS)", "1.0")
            rod_corr = self._prompt_float("Rod Length Correction (CR)", "1.0")
            gwt_depth = self._prompt_float("Groundwater table in feet", "2")
            unit_weight = self._prompt_float("Unit Weight of Soil (pcf)", "120")
            atm_pressure = self._prompt_float("Atmospheric Pressure (psf)", "2116")

            if hammer_eff <= 0 or hammer_eff > 1:
                print("Error: Hammer efficiency must be between 0 and 1")
                return
            if unit_weight <= 0:
                print("Error: Unit weight must be greater than 0")
                return
            if atm_pressure <= 0:
                print("Error: Atmospheric pressure must be greater than 0")
                return

            for boring_id, data in self.borings_data.items():
                print(f"\n{'=' * 110}")
                print(f"BORING: {boring_id}")
                print(f"{'=' * 110}")
                header = f"{'Depth':<8} {'N-fld':<8} {chr(963)+chr(39)+'v':<10} {'ER':<6} {'CB':<6} {'CS':<6} {'CR':<6} {'CN':<6} {'N60':<10} {'N160':<10}"
                print(header)
                print("-" * 110)

                for depth, n_field in zip(data['depths'], data['n_values']):
                    if gwt_depth >= depth:
                        eff_overburden = unit_weight * depth
                    else:
                        eff_overburden = (unit_weight * gwt_depth) + ((unit_weight - 62.4) * (depth - gwt_depth))
                    ER = hammer_eff / 0.6

                    n60 = convert_spt_to_n60(
                        n_field=n_field,
                        hammer_efficiency=hammer_eff,
                        borehole_diameter=borehole_corr,
                        sampler_correction=sampler_corr,
                        rod_length_correction=rod_corr
                    )

                    n160 = convert_spt_to_n160(
                        n_field=n_field,
                        hammer_efficiency=hammer_eff,
                        overburden_pressure=eff_overburden,
                        atmospheric_pressure=atm_pressure,
                        borehole_diameter=borehole_corr,
                        sampler_correction=sampler_corr,
                        rod_length_correction=rod_corr
                    )

                    if 'n60_values' not in data:
                        data['n60_values'] = []
                    if 'n160_values' not in data:
                        data['n160_values'] = []
                    data['n60_values'].append(n60)
                    data['n160_values'].append(n160)

                    CN = (atm_pressure / eff_overburden) ** 0.5

                    print(f"{depth:<8.1f} {n_field:<8} {eff_overburden:<10.2f} {ER:<6.3f} {borehole_corr:<6.2f} {sampler_corr:<6.2f} {rod_corr:<6.2f} {CN:<6.3f} {n60:<10.2f} {n160:<10.2f}")

            print(f"\n{'=' * 110}")
            print("CORRECTION PARAMETERS USED:")
            print(f"{'=' * 110}")
            print("Energy Correction:")
            print(f"  Hammer Efficiency: {hammer_eff:.2f} ({hammer_eff * 100:.0f}%)")
            print(f"  Energy Ratio (ER): {hammer_eff / 0.6:.3f}")
            print("\nEquipment Corrections:")
            print(f"  Borehole Diameter (CB): {borehole_corr:.2f}")
            print(f"  Sampler (CS): {sampler_corr:.2f}")
            print(f"  Rod Length (CR): {rod_corr:.2f}")
            print("\nOverburden Correction:")
            print(f"  Unit Weight: {unit_weight:.1f} pcf")
            print(f"  Atmospheric Pressure: {atm_pressure:.1f} psf")
            print("\nFormulas:")
            print("  N60 = N-field * ER * CB * CS * CR")
            print("  CN = (Pa / sigma'v)^0.5")
            print("  N160 = N60 * CN")
            print("  sigma'v = Unit Weight * Depth")
            print("\nNote: Depth units in feet, sigma'v in psf")

        except ValueError:
            print("Error: Invalid input values")

class SPT:

    def __init__(self):
        self.borings_data = {}

    def run(self):
        print("\n" + "=" * 70)
        print("Boring Information (Stratum Info)")
        print("=" * 70)

        while True:
            print("\nOptions:")
            print("  1. Add Boring")
            print("  2. View Soil Type Codes")
            print("  3. Clear All")
            print("  4. List Current Borings")
            print("  5. Proceed to N60 Calculator")
            print("  6. Back to Main Menu")
            choice = input("Enter choice: ").strip()

            if choice == '1':
                self.add_boring()
            elif choice == '2':
                self.soil_type_text_box()
            elif choice == '3':
                self.clear_all()
            elif choice == '4':
                self.list_borings()
            elif choice == '5':
                self.open_n60()
                return
            elif choice == '6':
                return
            else:
                print("Invalid choice")

    def open_n60(self):
        if not self.borings_data:
            print("Error: Please add at least one boring before proceeding to N60 calculator")
            return
        N60(self.borings_data).run()

    def _prompt(self, label, default):
        value = input(f"{label} [{default}]: ").strip()
        return value if value else default

    def list_borings(self):
        if not self.borings_data:
            print("No borings stored.")
            return
        print("\nStored borings:")
        for bid, data in self.borings_data.items():
            print(f"  {bid}: depths={data['depths']}, N={data['n_values']}, soil_types={data['soil types']}")

    def add_boring(self):
        boring_id = self._prompt("Boring ID", "B-1")
        n_values_str = self._prompt("N-values (comma-separated)", "10,29,40,17")
        soil_type_text = self._prompt("Soil type (comma-separated, 1-6)", "1,2,1,1")
        depths_str = self._prompt("Depths in feet (comma-separated)", "1,3,5,7")

        if not boring_id:
            print("Error: Please enter a Boring ID")
            return

        if not n_values_str or not depths_str:
            print("Error: Please enter both N-values and depths")
            return

        try:
            soil_types = [int(x.strip()) for x in soil_type_text.split(',')]
            for value in soil_types:
                if value < 1 or value > 6:
                    print("Error: Soil type must be between 1 and 6")
                    return
        except ValueError:
            print("Error: Invalid soil type value")
            return

        try:
            n_values = [int(x.strip()) for x in n_values_str.split(',')]
            depths = [float(x.strip()) for x in depths_str.split(',')]

            if len(n_values) != len(depths) or len(n_values) != len(soil_types):
                print("Error: Number of values must match")
                return

            self.borings_data[boring_id] = {
                'n_values': n_values,
                'depths': depths,
                'soil types': soil_types
            }

            print(f"Success: Boring {boring_id} added with {len(n_values)} readings")

        except ValueError:
            print("Error: Invalid input. Please use numbers only, separated by commas")

    def clear_all(self):
        self.borings_data = {}
        print("All data has been cleared")

    def soil_type_text_box(self):
        print("\nSoil Types:")
        print("  1: Sand (NC)")
        print("  2: Sand (saturated)")
        print("  3: Sand (OC)")
        print("  4: Gravelly Sand")
        print("  5: Clayey Sand")
        print("  6: Silts, sandy silts or clayey silt")

def main_menu():
    """Display the terminal main menu and dispatch to calculators."""
    while True:
        print("\n" + "=" * 70)
        print("SPT Calculator - Main Menu")
        print("=" * 70)
        print("\nSelect an option:")
        print("  1. Subsurface Information (SPT)")
        print("  2. N60 Calculator")
        print("  3. Drained Modulus Stiffness Calculator")
        print("  4. Footing Settlement on Sand")
        print("  5. Liquefaction Analysis (Youd et al. 2001)")
        print("  6. Exit")
        choice = input("Enter choice: ").strip()

        if choice == '1':
            SPT().run()
        elif choice == '2':
            N60().run()
        elif choice == '3':
            DrainedModulus().run()
        elif choice == '4':
            FootingSettlementSand({}).run()
        elif choice == '5':
            Liquefaction({}, magnitude=7.5, peak_acceleration_g=0.20, gwt_depth=5.0).run()
        elif choice == '6':
            print("Exiting...")
            return
        else:
            print("Invalid choice")

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
    # TODO: check these correlations and have the code return these in the list
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

class FootingSettlementSand:
    """
    Settlement of shallow spread and continuous footings bearing on sand.

    Two methods are evaluated per the 1996 revised Peck & Terzaghi practice
    manual:
      - Schmertmann (1978, revised 1996): strain-influence factor method
        with depth correction C1 and creep correction C2.
      - Burland and Burbridge (1985, revised 1996): empirical compressibility
        index Ic correlated to mean N60 in the depth of influence, with shape
        factor f_s, layer-thickness factor f_l, and time correction f_t.

    Soil data is consumed from the borings_data dict produced by SPT, N60,
    and DrainedModulus. Layer unit weights are derived through
    Stratums_and_SoilProps from each boring's N-values and soil types.

    Each method runs for every boring; worst/best/average borings are
    identified and differential settlement is reported as (max - min) and
    (max - avg). Settlement is reported at every user-specified time, with
    creep extension via Schmertmann's C2 and Burland-Burbridge's f_t.

    Footing shapes: square, rectangular, strip/continuous, circular.
    Load is supplied as gross applied pressure q (psf); effective overburden
    is subtracted internally using unit weights and the groundwater table.
    For Burland-Burbridge OC sand, supply the preconsolidation pressure
    sigma_p (psf); otherwise the sand is treated as NC.

    Can be invoked programmatically (full parameters via __init__) or
    interactively via run(), which prompts for footing geometry and load.
    """

    # Burland-Burbridge static-loading time-correction constants
    BB_R3 = 0.3
    BB_RT = 0.2

    # Unit conversions
    KPA_PER_PSF = 0.04788
    M_PER_FT = 0.3048

    def __init__(self,
                 borings_data: dict,
                 footing_width: Optional[float] = None,
                 footing_length: Optional[float] = None,
                 footing_shape: str = 'rectangular',
                 footing_depth: float = 0.0,
                 gross_pressure: float = 0.0,
                 gwt_depth: float = float('inf'),
                 sigma_p_psf: Optional[float] = None,
                 times_years: Tuple[float, ...] = (0.1, 1.0, 10.0, 50.0)):
        self.borings_data = borings_data
        self.B = footing_width
        self.L = footing_length if footing_length is not None else footing_width
        self.shape = (footing_shape or 'rectangular').lower()
        self.D = footing_depth
        self.q = gross_pressure
        self.gwt_depth = gwt_depth
        self.sigma_p_psf = sigma_p_psf
        self.times_years = tuple(times_years)
        self._normalize_shape()

        self.stratums = {}
        for bid, data in borings_data.items():
            self.stratums[bid] = self._build_stratums(data)

    def _normalize_shape(self):
        if self.B is None:
            return
        if self.shape in ('square', 'circular'):
            self.L = self.B
        elif self.shape in ('strip', 'continuous'):
            self.L = 1e6 * self.B   # plane-strain proxy: L/B >> 10

    def _build_stratums(self, data):
        depths = data.get('depths', [])
        n_for_corr = data.get('n_values') or data.get('n60_values') or []
        soil_types = data.get('soil types') or data.get('soil_types') or []
        if not (depths and n_for_corr and soil_types):
            return []
        layers = []
        for i, depth in enumerate(depths):
            thickness = depth if i == 0 else depths[i] - depths[i - 1]
            layers.append(Stratums_and_SoilProps(
                thickness=thickness,
                n_value=n_for_corr[i],
                soil_type=soil_types[i]
            ))
        return layers

    def _effective_overburden(self, boring_id: str, depth_ft: float) -> float:
        """Effective vertical stress (psf) at the given depth below ground surface."""
        if depth_ft <= 0:
            return 0.0
        layers = self.stratums.get(boring_id, [])
        if not layers:
            return 0.0
        sigma = 0.0
        cum = 0.0
        for layer in layers:
            top = cum
            bot = cum + layer.thickness
            slice_top = top
            slice_bot = min(depth_ft, bot)
            if slice_bot <= slice_top:
                break
            if slice_bot <= self.gwt_depth:
                sigma += layer.unit_weight * (slice_bot - slice_top)
            elif slice_top >= self.gwt_depth:
                sigma += (layer.unit_weight - 62.4) * (slice_bot - slice_top)
            else:
                above = self.gwt_depth - slice_top
                below = slice_bot - self.gwt_depth
                sigma += layer.unit_weight * above + (layer.unit_weight - 62.4) * below
            cum = bot
            if cum >= depth_ft:
                break
        return sigma

    def _schmertmann_iz_profile(self) -> Tuple[float, float, float]:
        """Return (Iz_top, z_peak_ft, z_zero_ft) interpolated by L/B."""
        ratio = self.L / self.B if self.B else 1.0
        if ratio < 1.0:
            ratio = 1.0
        # Schmertmann's linear interpolation: L/B=1 (axisym) -> L/B>=10 (plane strain)
        Iz_top = min(0.2, 0.1 + 0.0111 * (ratio - 1.0))
        z_peak = min(self.B, self.B * (0.5 + 0.0555 * (ratio - 1.0)))
        z_zero = min(4.0 * self.B, self.B * (2.0 + 0.222 * (ratio - 1.0)))
        return Iz_top, z_peak, z_zero

    def calculate_schmertmann(self, boring_id: str) -> dict:
        """Schmertmann settlement for one boring. Settlement reported in inches per time."""
        if boring_id not in self.borings_data:
            raise ValueError(f"Boring '{boring_id}' not found in borings_data")
        data = self.borings_data[boring_id]
        E_vals = data.get('E')
        if not E_vals:
            raise ValueError(f"Boring '{boring_id}' has no E values; run DrainedModulus first.")
        depths = data['depths']

        sigma_at_D = self._effective_overburden(boring_id, self.D)
        q_net = self.q - sigma_at_D

        Iz_top, z_peak, z_zero = self._schmertmann_iz_profile()

        if q_net <= 0:
            return {
                'boring_id': boring_id,
                'q_gross_psf': self.q,
                'sigma_v_at_D_psf': sigma_at_D,
                'q_net_psf': q_net,
                'note': 'q_net <= 0; no settlement.',
                'settlement_by_time_in': {t: 0.0 for t in self.times_years},
            }

        sigma_at_zpeak = self._effective_overburden(boring_id, self.D + z_peak)
        Iz_p = 0.5 + 0.1 * math.sqrt(q_net / max(sigma_at_zpeak, 1e-3))

        n = len(depths)
        bounds = []
        for i in range(n):
            top = 0.0 if i == 0 else 0.5 * (depths[i - 1] + depths[i])
            bot = (depths[i] + (depths[i] - top)) if i == n - 1 else 0.5 * (depths[i] + depths[i + 1])
            bounds.append((top, bot))

        sum_term = 0.0
        layer_breakdown = []
        for i in range(n):
            layer_top, layer_bot = bounds[i]
            z_t = layer_top - self.D
            z_b = layer_bot - self.D
            if z_b <= 0:
                continue
            if z_t >= z_zero:
                break
            z_t = max(0.0, z_t)
            z_b = min(z_zero, z_b)
            if z_b <= z_t:
                continue
            z_mid = 0.5 * (z_t + z_b)
            dz = z_b - z_t
            if z_mid <= z_peak:
                Iz = Iz_top + (Iz_p - Iz_top) * (z_mid / z_peak if z_peak > 0 else 0.0)
            else:
                Iz = Iz_p * (1.0 - (z_mid - z_peak) / (z_zero - z_peak)) if z_zero > z_peak else 0.0
            E_ksf = E_vals[i]
            if E_ksf is None or E_ksf <= 0:
                continue
            E_psf = E_ksf * 1000.0
            contribution = Iz * dz / E_psf
            sum_term += contribution
            layer_breakdown.append({
                'z_below_footing_ft': z_mid,
                'dz_ft': dz,
                'Iz': Iz,
                'E_ksf': E_ksf,
                'contribution': contribution,
            })

        C1 = max(0.5, 1.0 - 0.5 * sigma_at_D / q_net)
        settlement_by_time_in = {}
        for t in self.times_years:
            C2 = 1.0 + 0.2 * math.log10(t / 0.1)
            settlement_by_time_in[t] = C1 * C2 * q_net * sum_term * 12.0  # ft -> in

        return {
            'boring_id': boring_id,
            'q_gross_psf': self.q,
            'sigma_v_at_D_psf': sigma_at_D,
            'q_net_psf': q_net,
            'Iz_top': Iz_top,
            'Iz_peak': Iz_p,
            'z_peak_ft': z_peak,
            'z_zero_ft': z_zero,
            'C1': C1,
            'sum_Iz_dz_over_E': sum_term,
            'layer_breakdown': layer_breakdown,
            'settlement_by_time_in': settlement_by_time_in,
        }

    def calculate_burland_burbridge(self, boring_id: str) -> dict:
        """Burland-Burbridge settlement for one boring. Settlement reported in inches per time."""
        if boring_id not in self.borings_data:
            raise ValueError(f"Boring '{boring_id}' not found in borings_data")
        data = self.borings_data[boring_id]
        n60 = data.get('n60_values')
        if not n60:
            raise ValueError(f"Boring '{boring_id}' has no n60_values; run N60 first.")
        depths = data['depths']

        sigma_at_D = self._effective_overburden(boring_id, self.D)
        q_net_psf = self.q - sigma_at_D
        if q_net_psf <= 0:
            return {
                'boring_id': boring_id,
                'q_gross_psf': self.q,
                'sigma_v_at_D_psf': sigma_at_D,
                'q_net_psf': q_net_psf,
                'note': 'q_net <= 0; no settlement.',
                'settlement_by_time_in': {t: 0.0 for t in self.times_years},
            }

        q_net_kpa = q_net_psf * self.KPA_PER_PSF
        B_m = self.B * self.M_PER_FT
        zI_m = B_m ** 0.75
        zI_ft = zI_m / self.M_PER_FT

        samples = [n for d, n in zip(depths, n60) if self.D < d <= self.D + zI_ft]
        if not samples:
            below = [(d, n) for d, n in zip(depths, n60) if d >= self.D]
            if not below:
                raise ValueError(
                    f"No N60 data below footing depth {self.D} ft in boring '{boring_id}'"
                )
            samples = [below[0][1]]
        N_avg = sum(samples) / len(samples)

        sigma_p_kpa = (self.sigma_p_psf * self.KPA_PER_PSF) if self.sigma_p_psf is not None else None

        Ic_NC = 1.71 / (N_avg ** 1.4)
        Ic_OC = 0.57 / (N_avg ** 1.4)

        if sigma_p_kpa is None:
            Ic_reported = Ic_NC
            condition = 'NC'
            oc_split = False
        elif q_net_kpa <= sigma_p_kpa:
            Ic_reported = Ic_OC
            condition = 'OC, q_net <= sigma_p'
            oc_split = False
        else:
            Ic_reported = None
            condition = 'OC, q_net > sigma_p (split)'
            oc_split = True

        ratio = (self.L / self.B) if self.B else 1.0
        if ratio < 1.0:
            ratio = 1.0
        f_s = (1.25 * ratio / (ratio + 0.25)) ** 2

        depth_to_bottom = depths[-1] - self.D
        Hs_ft = max(0.0, min(depth_to_bottom, zI_ft))
        if Hs_ft >= zI_ft - 1e-9:
            f_l = 1.0
        else:
            f_l = (Hs_ft / zI_ft) * (2.0 - Hs_ft / zI_ft)

        settlement_by_time_in = {}
        for t in self.times_years:
            f_t = 1.0 + self.BB_R3 if t <= 3.0 else 1.0 + self.BB_R3 + self.BB_RT * math.log10(t / 3.0)
            if oc_split:
                S_mm = ((q_net_kpa - (2.0 / 3.0) * sigma_p_kpa) * Ic_NC
                        + (2.0 / 3.0) * sigma_p_kpa * Ic_OC) * (B_m ** 0.7) * f_s * f_l * f_t
            else:
                S_mm = q_net_kpa * (B_m ** 0.7) * Ic_reported * f_s * f_l * f_t
            settlement_by_time_in[t] = S_mm / 25.4

        return {
            'boring_id': boring_id,
            'q_gross_psf': self.q,
            'sigma_v_at_D_psf': sigma_at_D,
            'q_net_psf': q_net_psf,
            'q_net_kpa': q_net_kpa,
            'B_m': B_m,
            'zI_m': zI_m,
            'zI_ft': zI_ft,
            'N_avg': N_avg,
            'Ic': Ic_reported,
            'Ic_NC': Ic_NC,
            'Ic_OC': Ic_OC,
            'condition': condition,
            'f_s': f_s,
            'f_l': f_l,
            'Hs_ft': Hs_ft,
            'settlement_by_time_in': settlement_by_time_in,
        }

    def calculate_all(self) -> dict:
        out = {'schmertmann': {}, 'burland_burbridge': {}}
        for bid in self.borings_data:
            try:
                out['schmertmann'][bid] = self.calculate_schmertmann(bid)
            except (ValueError, KeyError, TypeError) as e:
                out['schmertmann'][bid] = {'error': str(e)}
            try:
                out['burland_burbridge'][bid] = self.calculate_burland_burbridge(bid)
            except (ValueError, KeyError, TypeError) as e:
                out['burland_burbridge'][bid] = {'error': str(e)}
        return out

    def calculate_differential(self) -> dict:
        all_results = self.calculate_all()
        differential = {}
        for method in ('schmertmann', 'burland_burbridge'):
            per_time = {}
            method_results = all_results[method]
            for t in self.times_years:
                settlements = {
                    bid: r['settlement_by_time_in'][t]
                    for bid, r in method_results.items()
                    if isinstance(r, dict) and 'settlement_by_time_in' in r
                }
                if not settlements:
                    continue
                worst_bid = max(settlements, key=settlements.get)
                best_bid = min(settlements, key=settlements.get)
                avg = sum(settlements.values()) / len(settlements)
                per_time[t] = {
                    'settlements_in': settlements,
                    'worst_boring': worst_bid,
                    'best_boring': best_bid,
                    'worst_settlement_in': settlements[worst_bid],
                    'best_settlement_in': settlements[best_bid],
                    'avg_settlement_in': avg,
                    'differential_max_minus_min_in': settlements[worst_bid] - settlements[best_bid],
                    'differential_max_minus_avg_in': settlements[worst_bid] - avg,
                }
            differential[method] = per_time
        return {'per_boring': all_results, 'differential': differential}

    def print_report(self) -> dict:
        results = self.calculate_differential()
        print("\n" + "=" * 92)
        print("FOOTING SETTLEMENT ON SAND")
        print("=" * 92)
        print(f"  Footing shape       : {self.shape}")
        print(f"  B = {self.B:.2f} ft, L = {self.L:.2f} ft, D = {self.D:.2f} ft")
        print(f"  Gross pressure q    : {self.q:.1f} psf")
        gwt_str = 'below influence zone' if self.gwt_depth == float('inf') else f"{self.gwt_depth:.1f} ft"
        print(f"  Groundwater table   : {gwt_str}")
        if self.sigma_p_psf is not None:
            print(f"  Preconsolidation    : sigma_p = {self.sigma_p_psf:.1f} psf (OC sand)")
        else:
            print(f"  Preconsolidation    : NC sand assumed")
        print(f"  Evaluation times    : {self.times_years} years")

        method_titles = {
            'schmertmann': "SCHMERTMANN (1978, revised 1996)",
            'burland_burbridge': "BURLAND-BURBRIDGE (1985, revised 1996)",
        }
        for method in ('schmertmann', 'burland_burbridge'):
            print("\n" + "-" * 92)
            print(method_titles[method])
            print("-" * 92)
            print(f"  {'Boring':<10}", end='')
            for t in self.times_years:
                print(f"  {('S('+str(t)+'yr) [in]'):<14}", end='')
            print()
            for bid, r in results['per_boring'][method].items():
                if 'error' in r:
                    print(f"  {bid:<10}  ERROR: {r['error']}")
                    continue
                print(f"  {bid:<10}", end='')
                for t in self.times_years:
                    val = r['settlement_by_time_in'].get(t, float('nan'))
                    print(f"  {val:<14.4f}", end='')
                print()

            diff = results['differential'].get(method, {})
            if diff:
                print(f"\n  Differential settlement across borings:")
                for t in self.times_years:
                    d = diff.get(t)
                    if not d:
                        continue
                    print(f"    t = {t} yr:")
                    print(f"      worst boring  : {d['worst_boring']:<6}  S = {d['worst_settlement_in']:.4f} in")
                    print(f"      best  boring  : {d['best_boring']:<6}  S = {d['best_settlement_in']:.4f} in")
                    print(f"      average across borings : {d['avg_settlement_in']:.4f} in")
                    print(f"      diff (max - min)       : {d['differential_max_minus_min_in']:.4f} in")
                    print(f"      diff (max - avg)       : {d['differential_max_minus_avg_in']:.4f} in")
        return results

    def run(self):
        """Interactive terminal flow. Expects borings_data already populated by SPT/N60/DrainedModulus."""
        if not self.borings_data:
            print("Error: no borings data. Run SPT -> N60 -> DrainedModulus first.")
            return
        missing = []
        for bid, data in self.borings_data.items():
            if not data.get('n60_values'):
                missing.append(f"{bid}: missing n60_values (run N60)")
            if not data.get('E'):
                missing.append(f"{bid}: missing E (run DrainedModulus)")
        if missing:
            print("Error: borings missing required data:")
            for m in missing:
                print(f"  {m}")
            return

        print("\n" + "=" * 70)
        print("Footing Settlement on Sand")
        print("=" * 70)
        try:
            shape_in = (input("Footing shape (1=square, 2=rectangular, 3=strip, 4=circular) [1]: ").strip()
                        or '1')
            shape_map = {'1': 'square', '2': 'rectangular', '3': 'strip', '4': 'circular'}
            if shape_in not in shape_map:
                print("Invalid shape.")
                return
            self.shape = shape_map[shape_in]

            self.B = float(input("Footing width B (ft) [5]: ").strip() or '5')
            if self.shape == 'rectangular':
                self.L = float(input(f"Footing length L (ft) [{self.B}]: ").strip() or str(self.B))
            else:
                self._normalize_shape()

            self.D = float(input("Footing depth D below ground (ft) [2]: ").strip() or '2')
            self.q = float(input("Gross applied pressure q (psf) [3000]: ").strip() or '3000')

            gwt_in = input("Groundwater table depth (ft) [blank for none]: ").strip()
            self.gwt_depth = float(gwt_in) if gwt_in else float('inf')

            oc_in = (input("Sand condition (n=NC, o=OC) [n]: ").strip().lower() or 'n')
            if oc_in.startswith('o'):
                self.sigma_p_psf = float(input("Preconsolidation pressure sigma_p (psf): ").strip())
            else:
                self.sigma_p_psf = None

            times_in = input("Evaluation times in years (comma-separated) [0.1,1,10,50]: ").strip() or "0.1,1,10,50"
            self.times_years = tuple(float(x.strip()) for x in times_in.split(','))
        except ValueError as e:
            print(f"Error: invalid input ({e})")
            return

        self.print_report()

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

    def calculate_rankine_ka_profile(self, boring_id: str, boring_data: dict, gamma: float, phi: float,
                                     alpha: float = 0) -> dict:
        """
        Calculate K_a (active earth pressure coefficient) for all depths in a boring
        for c-phi soil with vertical backface

        Parameters:
        boring_id: ID of the boring to analyze
        boring_data: Dictionary containing boring information
        gamma: unit weight of soil
        phi: internal friction angle (in radians)
        alpha: slope angle (in radians)

        Returns:
        Dictionary with depths and corresponding K_a values
        """
        if boring_id not in boring_data:
            raise ValueError(f"Boring ID {boring_id} not found in boring_data")

        boring = boring_data[boring_id]
        depths = boring['depths']
        cohesion_values = boring['cohesion']

        ka_values = []
        z_crit_values = []

        for i, depth in enumerate(depths):
            c = cohesion_values[i]
            gamma_z = gamma * depth

            # Handle case where depth is 0 or very small
            if gamma_z < 1e-6:
                ka_values.append(None)  # or use a default value
                continue

            cos_phi = np.cos(phi)
            cos_alpha = np.cos(alpha)
            sin_phi = np.sin(phi)

            # Term inside the square root
            sqrt_term = np.sqrt(
                4 * cos_alpha ** 2 * (cos_alpha ** 2 - cos_phi ** 2) +
                4 * (c / gamma_z) ** 2 * cos_phi ** 2 +
                8 * (c / gamma_z) * cos_alpha ** 2 * sin_phi * cos_phi
            )

            # Main calculation
            K_a = (1 / cos_phi ** 2) * (
                    2 * cos_alpha ** 2 + 2 * (c / gamma_z) * cos_phi * sin_phi -
                    sqrt_term - 1
            )

            z_crit = ((2 * c) / gamma) * (np.sqrt((1 + np.sin(phi)) / (1 - np.sin(phi))))

            z_crit_values.append(z_crit)
            ka_values.append(K_a)

        return {
            'boring_id': boring_id,
            'depths': depths,
            'ka_values': ka_values,
            'cohesion': cohesion_values,
            'critical depth': z_crit_values
        }

    def calculate_rankine_ka_all_borings(self, boring_data: dict, gamma: float, phi: float, alpha: float = 0) -> dict:
        """
        Calculate K_a for all borings in the dataset

        Parameters:
        boring_data: Dictionary containing all boring information
        gamma: unit weight of soil
        phi: internal friction angle (in radians)
        alpha: slope angle (in radians)

        Returns:
        Dictionary with boring_id as keys and Ka profiles as values
        """
        results = {}

        for boring_id in boring_data.keys():
            results[boring_id] = self.calculate_rankine_ka_profile(
                boring_id, boring_data, gamma, phi, alpha
            )

        return results
    
    def calculate_rankine_kp_profile(self, boring_id: str, boring_data: dict, gamma: float, phi: float,
                                     alpha: float = 0) -> dict:
        """
        Calculate K_p (passive earth pressure coefficient) for all depths in a boring
        for c-phi soil with vertical backface

        Parameters:
        boring_id: ID of the boring to analyze
        boring_data: Dictionary containing boring information
        gamma: unit weight of soil
        phi: internal friction angle (in radians)
        alpha: slope angle (in radians)

        Returns:
        Dictionary with depths and corresponding K_a values
        """
        if boring_id not in boring_data:
            raise ValueError(f"Boring ID {boring_id} not found in boring_data")

        boring = boring_data[boring_id]
        depths = boring['depths']
        cohesion_values = boring['cohesion']

        ka_values = []
        z_crit_values = []

        for i, depth in enumerate(depths):
            c = cohesion_values[i]
            gamma_z = gamma * depth

            # Handle case where depth is 0 or very small
            if gamma_z < 1e-6:
                ka_values.append(None)  # or use a default value
                continue

            cos_phi = np.cos(phi)
            cos_alpha = np.cos(alpha)
            sin_phi = np.sin(phi)

            # Term inside the square root
            sqrt_term = np.sqrt(
                4 * cos_alpha ** 2 * (cos_alpha ** 2 - cos_phi ** 2) +
                4 * (c / gamma_z) ** 2 * cos_phi ** 2 +
                8 * (c / gamma_z) * cos_alpha ** 2 * sin_phi * cos_phi
            )

            K_a = (1 / cos_phi ** 2) * (
                    2 * cos_alpha ** 2 + 2 * (c / gamma_z) * cos_phi * sin_phi -
                    sqrt_term - 1
            )
            K_p = 1 / K_a if K_a != 0 else None

            z_crit = ((2 * c) / gamma) * (np.sqrt((1 + np.sin(phi)) / (1 - np.sin(phi))))

            z_crit_values.append(z_crit)
            ka_values.append(K_p)

        return {
            'boring_id': boring_id,
            'depths': depths,
            'ka_values': ka_values,
            'cohesion': cohesion_values,
            'critical depth': z_crit_values
        }
    def calculate_rankine_kp_all_borings(self, boring_data: dict, gamma: float, phi: float, alpha: float = 0) -> dict:
        """
        Calculate K_p for all borings in the dataset

        Parameters:
        boring_data: Dictionary containing all boring information
        gamma: unit weight of soil
        phi: internal friction angle (in radians)
        alpha: slope angle (in radians)

        Returns:
        Dictionary with boring_id as keys and Kp profiles as values
        """
        results = {}

        for boring_id in boring_data.keys():
            results[boring_id] = self.calculate_rankine_kp_profile(
                boring_id, boring_data, gamma, phi, alpha
            )

        return results


    @staticmethod
    def _interpolate_critical_wedge(trials: List[dict], force_key: str,
                                    kind: str) -> Tuple[float, float, dict, List[str]]:
        """
        Find the critical trial wedge by quadratic interpolation across trials.

        kind='max' selects the wedge with the largest force (active case);
        kind='min' selects the smallest (passive case). Trials with a near-zero
        denominator or non-physical force (negative or non-finite) are dropped
        before fitting. If the interpolated vertex falls outside the surviving
        trial theta range, or the parabola has the wrong concavity, the result
        falls back to the trial extreme.

        Returns (theta_crit_deg, force_crit, quadratic_fit_dict, warnings).
        """
        assert kind in ('max', 'min')
        warnings: List[str] = []

        valid = [t for t in trials
                 if math.isfinite(t[force_key]) and t[force_key] > 0
                 and t[force_key] < 1e12]
        if len(valid) < len(trials):
            warnings.append(
                f"Dropped {len(trials) - len(valid)} trial(s) with non-physical "
                f"{force_key} (singular denominator or out of range)."
            )

        thetas_valid = np.array([t['theta_deg'] for t in valid])
        forces_valid = np.array([t[force_key] for t in valid])

        extreme_idx = (int(np.argmax(forces_valid)) if kind == 'max'
                       else int(np.argmin(forces_valid)))
        extreme_theta = float(thetas_valid[extreme_idx])
        extreme_force = float(forces_valid[extreme_idx])

        if len(valid) < 3:
            warnings.append("Fewer than 3 valid trials — skipping quadratic fit; "
                            "reporting the extreme trial directly.")
            fit = {'a': float('nan'), 'b': float('nan'), 'c': float('nan')}
            return extreme_theta, extreme_force, fit, warnings

        a, b, c_const = np.polyfit(thetas_valid, forces_valid, 2)
        fit = {'a': float(a), 'b': float(b), 'c': float(c_const)}

        concavity_ok = (a < 0) if kind == 'max' else (a > 0)
        if not concavity_ok:
            warnings.append(
                "Quadratic concavity does not bracket a true "
                f"{'maximum' if kind == 'max' else 'minimum'}; "
                "falling back to the extreme trial."
            )
            return extreme_theta, extreme_force, fit, warnings

        theta_crit = float(-b / (2.0 * a))
        if not (thetas_valid.min() <= theta_crit <= thetas_valid.max()):
            warnings.append(
                f"Interpolated vertex theta={theta_crit:.2f} deg falls outside the "
                f"valid trial range [{thetas_valid.min():.2f}, {thetas_valid.max():.2f}]; "
                "falling back to the extreme trial."
            )
            return extreme_theta, extreme_force, fit, warnings

        force_crit = float(a * theta_crit ** 2 + b * theta_crit + c_const)
        return theta_crit, force_crit, fit, warnings

    def _wedge_geometry(self, H: float, theta_rad: float,
                        alpha_rad: float, beta_rad: float) -> Tuple[float, float, float, float]:
        """
        Compute trial-wedge geometry for a battered wall with sloped backfill.

        Coordinate frame: wall heel at origin; +x into backfill; +y up.
        Wall back runs from (0, 0) to (-H * tan alpha, H), so positive alpha
        leans the wall AWAY from the backfill at the top. Backfill rises from
        the top of the wall at angle beta from horizontal. The failure plane
        leaves the heel at angle theta from horizontal.

        Returns (x_top, y_top, area, L_failure) where (x_top, y_top) is the
        intersection of the failure plane with the backfill surface.
        """
        tan_a = math.tan(alpha_rad)
        tan_b = math.tan(beta_rad)
        tan_t = math.tan(theta_rad)

        denom = tan_t - tan_b
        if denom <= 0:
            raise ValueError(
                f"Failure plane angle theta={math.degrees(theta_rad):.2f} deg must exceed "
                f"backfill slope beta={math.degrees(beta_rad):.2f} deg"
            )

        x_top = H * (1.0 + tan_a * tan_b) / denom
        y_top = x_top * tan_t

        # Shoelace area for vertices (0,0), (-H*tan_a, H), (x_top, y_top).
        area = 0.5 * abs(H * (x_top + tan_a * y_top))
        L_failure = math.hypot(x_top, y_top)

        return x_top, y_top, area, L_failure

    def calculate_coulomb_active_wedge(self, wall_height: float, gamma: float,
                                       phi: float, c: float = 0.0,
                                       delta: Optional[float] = None,
                                       wall_batter_deg: float = 0.0,
                                       backfill_slope_deg: float = 0.0) -> dict:
        """
        Compute the active earth thrust using a trial-wedge force-polygon method.

        Geometry (see _wedge_geometry for the full coordinate convention):
          - Vertical retained height H.
          - wall_batter_deg (alpha): wall-back angle from vertical;
              positive => wall leans away from backfill at the top.
          - backfill_slope_deg (beta): backfill angle from horizontal;
              positive => rising backfill.

        Five trial failure surfaces are rotated through the soil starting from
        the Rankine angle theta0 = 45 + phi/2 (measured CCW from horizontal at
        the heel) at offsets -10, -5, 0, +5, +10 degrees.

        Forces on the wedge (per linear foot of wall):
            W  : self-weight, gamma * A                                [down]
            C  : cohesion along failure plane, c * L                   [up-slope]
            R  : reaction on failure plane, inclined phi from the failure-plane
                 normal on the up-slope side (direction known, magnitude unknown)
            Pa : reaction from wall on wedge, inclined delta above the wall's
                 outward normal (so at angle alpha + delta from horizontal)
                 (direction known, magnitude unknown)

        Closing the polygon (x- and y-equilibrium) yields, in closed form:
            Pa = [W * sin(theta - phi) - C * cos(phi)]
                 / cos(theta - phi - delta - alpha)

        As theta varies, Pa(theta) traces a curve whose peak is the critical
        active wedge (largest force the wall must resist). A quadratic is fit
        through the 5 trial points and the vertex (dPa/dtheta = 0) is taken as
        the interpolated critical wedge.

        Args:
            wall_height:        H, vertical retained height
            gamma:              soil unit weight
            phi:                internal friction angle (degrees)
            c:                  cohesion on failure plane
            delta:              wall friction angle (degrees); default (2/3)*phi
            wall_batter_deg:    alpha (degrees), default 0 (vertical wall)
            backfill_slope_deg: beta (degrees), default 0 (horizontal backfill)

        Returns:
            Dict with each trial's geometry/forces, the initial Rankine angle,
            the interpolated critical wedge angle and Pa, plus the input geometry
            for traceability.
        """
        H = wall_height
        phi_rad = math.radians(phi)
        if delta is None:
            delta = (2.0 / 3.0) * phi
        delta_rad = math.radians(delta)
        alpha_rad = math.radians(wall_batter_deg)
        beta_rad = math.radians(backfill_slope_deg)

        theta0_rad = math.radians(45.0) + phi_rad / 2.0

        trials = []
        for offset_deg in (-10, -5, 0, 5, 10):
            theta_rad = theta0_rad + math.radians(offset_deg)
            x_top, y_top, area, L = self._wedge_geometry(H, theta_rad, alpha_rad, beta_rad)

            W = gamma * area
            Cf = c * L

            num = W * math.sin(theta_rad - phi_rad) - Cf * math.cos(phi_rad)
            den = math.cos(theta_rad - phi_rad - delta_rad - alpha_rad)
            Pa = num / den if abs(den) > 1e-9 else float('inf')

            sin_tp = math.sin(theta_rad - phi_rad)
            R = ((Cf * math.cos(theta_rad) + Pa * math.cos(alpha_rad + delta_rad)) / sin_tp
                 if sin_tp != 0 and math.isfinite(Pa) else float('inf'))

            trials.append({
                'offset_deg': offset_deg,
                'theta_deg': math.degrees(theta_rad),
                'L_failure': L,
                'wedge_area': area,
                'wedge_weight': W,
                'cohesion_force': Cf,
                'R': R,
                'Pa': Pa,
            })

        theta_crit_deg, pa_crit, fit, warns = self._interpolate_critical_wedge(
            trials, force_key='Pa', kind='max')

        return {
            'trials': trials,
            'theta0_deg': math.degrees(theta0_rad),
            'theta_critical_deg': theta_crit_deg,
            'Pa_critical': pa_crit,
            'delta_deg': delta,
            'wall_batter_deg': wall_batter_deg,
            'backfill_slope_deg': backfill_slope_deg,
            'quadratic_fit': fit,
            'warnings': warns,
        }

    def calculate_coulomb_passive_wedge(self, wall_height: float, gamma: float,
                                        phi: float, c: float = 0.0,
                                        delta: Optional[float] = None,
                                        wall_batter_deg: float = 0.0,
                                        backfill_slope_deg: float = 0.0) -> dict:
        """
        Compute the passive earth thrust using a trial-wedge force-polygon method.

        Geometry conventions match calculate_coulomb_active_wedge.

        For passive failure the wedge slides UP the failure plane (away from the
        wall), so friction on the failure plane and on the wall both reverse
        direction relative to the active case. The reaction R is now inclined phi
        from the failure-plane normal on the DOWN-slope side, and Pp is inclined
        delta BELOW the wall's outward normal (at angle alpha - delta from
        horizontal). Cohesion still opposes wedge motion, so it now points
        down-slope (and adds to the resisting force).

        Trial failure surfaces are rotated through the soil starting from the
        Rankine passive angle theta0 = 45 - phi/2, at offsets -10, -5, 0, +5, +10
        degrees. The polygon closes for:
            Pp = [W * sin(theta + phi) + C * cos(phi)]
                 / cos(theta + phi + delta - alpha)

        The critical passive wedge MINIMISES Pp (smallest force the wall can
        mobilise to resist motion). A quadratic is fit through the 5 trial
        points; the vertex of the upward-opening parabola is the interpolated
        critical wedge.

        Returns:
            Dict with each trial's geometry/forces, theta0, theta_critical, the
            critical Pp, and the input geometry for traceability.
        """
        H = wall_height
        phi_rad = math.radians(phi)
        if delta is None:
            delta = (2.0 / 3.0) * phi
        delta_rad = math.radians(delta)
        alpha_rad = math.radians(wall_batter_deg)
        beta_rad = math.radians(backfill_slope_deg)

        theta0_rad = math.radians(45.0) - phi_rad / 2.0

        trials = []
        for offset_deg in (-10, -5, 0, 5, 10):
            theta_rad = theta0_rad + math.radians(offset_deg)
            x_top, y_top, area, L = self._wedge_geometry(H, theta_rad, alpha_rad, beta_rad)

            W = gamma * area
            Cf = c * L

            num = W * math.sin(theta_rad + phi_rad) + Cf * math.cos(phi_rad)
            den = math.cos(theta_rad + phi_rad + delta_rad - alpha_rad)
            Pp = num / den if abs(den) > 1e-9 else float('inf')

            sin_tp = math.sin(theta_rad + phi_rad)
            R = ((Pp * math.cos(alpha_rad - delta_rad) - Cf * math.cos(theta_rad)) / sin_tp
                 if sin_tp != 0 and math.isfinite(Pp) else float('inf'))

            trials.append({
                'offset_deg': offset_deg,
                'theta_deg': math.degrees(theta_rad),
                'L_failure': L,
                'wedge_area': area,
                'wedge_weight': W,
                'cohesion_force': Cf,
                'R': R,
                'Pp': Pp,
            })

        theta_crit_deg, pp_crit, fit, warns = self._interpolate_critical_wedge(
            trials, force_key='Pp', kind='min')

        return {
            'trials': trials,
            'theta0_deg': math.degrees(theta0_rad),
            'theta_critical_deg': theta_crit_deg,
            'Pp_critical': pp_crit,
            'delta_deg': delta,
            'wall_batter_deg': wall_batter_deg,
            'backfill_slope_deg': backfill_slope_deg,
            'quadratic_fit': fit,
            'warnings': warns,
        }

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

                        # Coulomb is now computed as a total-thrust trial-wedge analysis
                        # via calculate_coulomb_active_wedge(); the per-layer K profile
                        # uses Rankine.
                        if side == 'active':
                            K = self.calculate_rankine_ka(layer.phi)
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

class LogSpiralEarthPressure:

    """
    Passive earth pressure via a log-spiral failure surface.

    Geometry (passive case, vertical wall, horizontal ground):
      - Top of wall at A=(0, 0) (ground level); bottom at B=(0, -H).
      - Soil in +x. The Rankine line through A is rotated (45 - phi/2)
        clockwise from horizontal — into the soil it has direction
        (cos α, -sin α), α = 45° - phi/2. The trial spiral centers are
        placed on the OPPOSITE side of A (above ground, behind the wall),
        i.e. at parameter d along direction (-cos α, sin α):
            O = (-d cos α, d sin α),  d >= 0.
      - Trials start at d = 0.10*H and step outward by 0.05*H per trial,
        moving O further away from the wall along that line, out to a
        maximum of d = 3.0*H (default num_trials = 59).
      - The spiral passes through B and is extended CCW until it intersects
        the Rankine line at R, on the soil side of O. The wedge is closed by
        the vertical line from R up to the ground at D=(x_R, 0), then back
        along the ground to A.

    Forces on the wedge (per unit wall length):
      Pp : passive thrust on the wall face, inclined δ BELOW horizontal,
           resultant applied H/3 above the wall bottom (triangular soil
           pressure assumption).
      W  : self-weight γ * Area, applied at the wedge centroid.
      PR : Rankine ACTIVE thrust on the vertical RD line (height h = |y_R|):
             PR = 0.5 γ h² Ka + q h Ka - 2 c h sqrt(Ka),
             Ka = (1-sinφ)/(1+sinφ);
           horizontal, acting on the wedge toward the wall, at the centroid
           of the active pressure distribution.
      C  : cohesion along the spiral surface (closed-form moment about O:
           M_c = -c (r_R² - r_B²) / (2 tan φ); degenerates to -c r² Δθ for
           φ → 0).
      Q  : surcharge q on the ground surface from x=0 to x=x_R, at midspan.
      Normal + Mohr-Coulomb friction on the spiral surface have a resultant
      passing through O (log-spiral property), contributing zero moment.

    Moment equilibrium about O (CCW positive) solves Pp directly. The
    critical wedge MINIMISES Pp across trials; a quadratic fit interpolates
    the minimum (falling back to the smallest trial Pp if the fit isn't
    upward-opening or the vertex falls outside the trial range).
    """


    def __init__(self,
                 wall_height: float,
                 gamma: float,
                 phi: float,
                 cohesion: float = 0.0,
                 wall_friction_delta: Optional[float] = None,
                 surcharge: float = 0.0,
                 num_trials: int = 59,
                 d_init_frac: float = 0.10,
                 d_step_frac: float = 0.05,
                 spiral_segments: int = 60):
        """
        Args:
            wall_height: H (ft or m)
            gamma: soil unit weight (pcf or kN/m³)
            phi: internal friction angle (degrees)
            cohesion: c (psf or kPa)
            wall_friction_delta: δ (degrees); default = (2/3) phi
            surcharge: q on ground surface (psf or kPa)
            num_trials: number of trial wedges
            d_init_frac: first trial d as fraction of H (default 0.10).
                d is measured along the Rankine line in the direction
                opposite the soil (above ground, behind the wall).
            d_step_frac: increment between trials as fraction of H (default 0.05).
                Each step moves O further away from the wall along that line.
            spiral_segments: spiral discretisation for area/centroid
        """
        self.H = wall_height
        self.gamma = gamma
        self.phi = phi
        self.c = cohesion
        self.delta = wall_friction_delta if wall_friction_delta is not None else (2.0/3.0) * phi
        self.q = surcharge
        self.num_trials = num_trials
        self.d_init = d_init_frac * wall_height
        self.d_step = d_step_frac * wall_height
        self.spiral_segments = spiral_segments

        self.phi_rad = math.radians(phi)
        self.delta_rad = math.radians(self.delta)
        self.alpha_rad = math.radians(45.0) - self.phi_rad / 2.0
        self.Ka = (1.0 - math.sin(self.phi_rad)) / (1.0 + math.sin(self.phi_rad))

    def _passive_trial(self, d: float) -> dict:
        """Compute Pp for one trial wedge with spiral center O at parameter d on Rankine line."""
        H = self.H
        phi_rad = self.phi_rad
        delta_rad = self.delta_rad
        alpha = self.alpha_rad
        gamma = self.gamma
        c = self.c
        q = self.q

        cos_a = math.cos(alpha)
        sin_a = math.sin(alpha)

        # Spiral center on Rankine line, OPPOSITE direction from the soil
        # (above ground, behind the wall). Each trial steps O further away.
        Ox = -d * cos_a
        Oy = d * sin_a

        # Radius and angle from O to wall bottom B = (0, -H)
        dx0 = -Ox
        dy0 = -H - Oy
        r_start = math.hypot(dx0, dy0)
        theta_start = math.atan2(dy0, dx0)

        # Rankine line crosses A and continues into the soil on the far side
        # of O at polar angle theta = -alpha (from O). Reach it CCW from B.
        theta_R = -alpha
        delta_theta = theta_R - theta_start
        if delta_theta <= 0:
            delta_theta += 2.0 * math.pi

        # Intersection point R lies along the Rankine direction at parameter
        # t_R = r_R - d from A (since O is at -d along the Rankine line).

        tan_phi = math.tan(phi_rad)
        if abs(tan_phi) < 1e-9:
            r_R = r_start
        else:
            r_R = r_start * math.exp(delta_theta * tan_phi)

        # With O at -d along the Rankine line, R = O + r_R*(cos α, -sin α),
        # which lies at parameter t_R = r_R - d from A along the soil side.
        t_R = r_R - d
        x_R = t_R * cos_a
        y_R = -t_R * sin_a
        h = -y_R   # height of Rankine zone (y_R is negative below ground)

        # Discretise spiral from B to R; build closed wedge polygon CCW:
        #   A=(0,0) -> B=(0,-H) -> spiral arc -> R=(x_R,y_R) -> D=(x_R,0) -> A
        N = self.spiral_segments
        spiral_pts = []
        for i in range(N + 1):
            th = theta_start + delta_theta * (i / N)
            r_i = r_start if abs(tan_phi) < 1e-9 else r_start * math.exp((th - theta_start) * tan_phi)
            spiral_pts.append((Ox + r_i * math.cos(th), Oy + r_i * math.sin(th)))

        poly = [(0.0, 0.0)]
        poly.extend(spiral_pts)   # spiral_pts[0] ≈ B, spiral_pts[-1] = R
        poly.append((x_R, 0.0))   # D

        # Shoelace area + centroid
        cross_sum = 0.0
        cx_num = 0.0
        cy_num = 0.0
        n_poly = len(poly)
        for i in range(n_poly):
            x1, y1 = poly[i]
            x2, y2 = poly[(i + 1) % n_poly]
            cross = x1 * y2 - x2 * y1
            cross_sum += cross
            cx_num += (x1 + x2) * cross
            cy_num += (y1 + y2) * cross
        area_signed = cross_sum / 2.0
        if abs(area_signed) < 1e-12:
            raise ValueError(f"Degenerate wedge at d={d}")
        Cx = cx_num / (6.0 * area_signed)
        Cy = cy_num / (6.0 * area_signed)
        area = abs(area_signed)

        W = gamma * area

        # Rankine active force on vertical line RD (height h)
        Ka = self.Ka
        sqrt_Ka = math.sqrt(Ka)
        PR = 0.5 * gamma * h * h * Ka + q * h * Ka - 2.0 * c * h * sqrt_Ka

        # Centroid of active pressure distribution (depth below ground = positive)
        M_about_top = (0.5 * q * Ka * h * h
                       + gamma * Ka * h * h * h / 3.0
                       - c * sqrt_Ka * h * h)
        depth_PR = M_about_top / PR if abs(PR) > 1e-9 else h / 2.0
        y_PR = -depth_PR

        # Moments about O, CCW positive
        # Pp at (0, -2H/3), F = (Pp cos δ, -Pp sin δ)
        Pp_arm = Ox * math.sin(delta_rad) + (2.0 * H / 3.0 + Oy) * math.cos(delta_rad)
        M_W = -W * (Cx - Ox)
        M_PR = PR * (y_PR - Oy)
        if abs(tan_phi) < 1e-9:
            M_c = -c * r_start * r_start * delta_theta
        else:
            M_c = -c * (r_R * r_R - r_start * r_start) / (2.0 * tan_phi)
        M_q = -q * x_R * (x_R / 2.0 - Ox)

        if abs(Pp_arm) < 1e-9:
            Pp = float('inf')
        else:
            Pp = -(M_W + M_PR + M_c + M_q) / Pp_arm

        return {
            'd': d,
            'O': (Ox, Oy),
            'r_B': r_start,
            'r_R': r_R,
            'theta_B_deg': math.degrees(theta_start),
            'theta_R_deg': math.degrees(theta_R),
            'spiral_sweep_deg': math.degrees(delta_theta),
            'intersection': (x_R, y_R),
            'rankine_height': h,
            'wedge_area': area,
            'centroid': (Cx, Cy),
            'W': W,
            'PR': PR,
            'PR_depth': depth_PR,
            'M_W': M_W,
            'M_PR': M_PR,
            'M_c': M_c,
            'M_q': M_q,
            'Pp_arm': Pp_arm,
            'Pp': Pp,
        }

    def calculate_critical_passive(self) -> dict:
        """Run trial wedges and interpolate the critical (minimum-Pp) wedge."""
        trials = []
        for i in range(self.num_trials):
            d = self.d_init + i * self.d_step
            try:
                trials.append(self._passive_trial(d))
            except (ValueError, ZeroDivisionError) as e:
                trials.append({'d': d, 'error': str(e), 'Pp': float('inf')})

        valid = [t for t in trials
                 if 'error' not in t and math.isfinite(t['Pp']) and t['Pp'] > 0]
        warnings: List[str] = []
        if not valid:
            return {
                'trials': trials,
                'critical_d': None,
                'Pp_critical': None,
                'quadratic_fit': None,
                'warnings': ['No valid trials produced a finite, positive Pp.'],
            }

        d_arr = np.array([t['d'] for t in valid])
        Pp_arr = np.array([t['Pp'] for t in valid])
        idx = int(np.argmin(Pp_arr))
        d_min = float(d_arr[idx])
        Pp_min = float(Pp_arr[idx])

        if len(valid) < 3:
            warnings.append("Fewer than 3 valid trials; reporting the minimum trial Pp.")
            return {'trials': trials, 'critical_d': d_min, 'Pp_critical': Pp_min,
                    'quadratic_fit': None, 'warnings': warnings}

        a, b, c_fit = np.polyfit(d_arr, Pp_arr, 2)
        fit = {'a': float(a), 'b': float(b), 'c': float(c_fit)}

        if a <= 0:
            warnings.append("Quadratic fit is not upward-opening; falling back to minimum trial.")
            return {'trials': trials, 'critical_d': d_min, 'Pp_critical': Pp_min,
                    'quadratic_fit': fit, 'warnings': warnings}

        d_crit = -b / (2.0 * a)
        if not (d_arr.min() <= d_crit <= d_arr.max()):
            warnings.append(
                f"Interpolated vertex d={d_crit:.3f} outside trial range "
                f"[{d_arr.min():.3f}, {d_arr.max():.3f}]; using minimum trial.")
            return {'trials': trials, 'critical_d': d_min, 'Pp_critical': Pp_min,
                    'quadratic_fit': fit, 'warnings': warnings}

        Pp_crit = float(a * d_crit * d_crit + b * d_crit + c_fit)
        if Pp_crit <= 0 or Pp_crit < 0.5 * Pp_min:
            warnings.append(
                f"Interpolated vertex Pp={Pp_crit:.2f} is below the data envelope "
                f"(min trial Pp={Pp_min:.2f}); trial range likely doesn't bracket a "
                "minimum (Pp(d) appears monotonic). Falling back to minimum trial.")
            return {'trials': trials, 'critical_d': d_min, 'Pp_critical': Pp_min,
                    'quadratic_fit': fit, 'warnings': warnings}
        return {'trials': trials, 'critical_d': float(d_crit),
                'Pp_critical': Pp_crit, 'quadratic_fit': fit, 'warnings': warnings}

    def print_report(self) -> dict:
        result = self.calculate_critical_passive()
        print("\n" + "=" * 110)
        print("LOG-SPIRAL PASSIVE EARTH PRESSURE")
        print("=" * 110)
        print(f"  Wall height H            : {self.H:.3f}")
        print(f"  Soil unit weight gamma   : {self.gamma:.2f}")
        print(f"  Friction angle phi       : {self.phi:.2f} deg")
        print(f"  Cohesion c               : {self.c:.2f}")
        print(f"  Wall friction delta      : {self.delta:.2f} deg")
        print(f"  Surcharge q              : {self.q:.2f}")
        print(f"  Rankine angle 45 - phi/2 : {math.degrees(self.alpha_rad):.2f} deg")
        print(f"  Ka                       : {self.Ka:.4f}")
        print(f"  Trials                   : {self.num_trials}, d_start={self.d_init:.3f}, "
              f"step={self.d_step:.3f}")
        print("-" * 110)
        print(f"  {'d':<8} {'Ox':<8} {'Oy':<9} {'x_R':<8} {'y_R':<9} {'h':<8} "
              f"{'Area':<10} {'W':<10} {'PR':<10} {'Pp':<12}")
        for t in result['trials']:
            if 'error' in t:
                print(f"  {t['d']:<8.3f} ERROR: {t['error']}")
                continue
            Ox, Oy = t['O']
            xR, yR = t['intersection']
            print(f"  {t['d']:<8.3f} {Ox:<8.3f} {Oy:<9.3f} {xR:<8.3f} {yR:<9.3f} "
                  f"{t['rankine_height']:<8.3f} {t['wedge_area']:<10.2f} "
                  f"{t['W']:<10.2f} {t['PR']:<10.2f} {t['Pp']:<12.2f}")
        print("-" * 110)
        if result.get('Pp_critical') is not None:
            print(f"  Critical wedge: d = {result['critical_d']:.4f}")
            print(f"  Critical passive force Pp = {result['Pp_critical']:.2f} per unit wall length")
        for w in result.get('warnings', []):
            print(f"  WARNING: {w}")
        print("=" * 110)
        return result


class Liquefaction:
    """
    Liquefaction triggering analysis and factor of safety per Youd et al. (2001),
    "Liquefaction Resistance of Soils: Summary Report from the 1996 NCEER and
    1998 NCEER/NSF Workshops on Evaluation of Liquefaction Resistance of Soils,"
    J. Geotech. Geoenviron. Eng., 127(10), 817-833.

    Procedure (SPT-based, simplified Seed-Idriss method):
      1) Cyclic Stress Ratio (CSR), Seed and Idriss (1971):
             CSR = 0.65 * (a_max/g) * (sigma_v / sigma'_v) * r_d
         Stress reduction r_d from Liao and Whitman (1986), as adopted in Youd
         et al. (2001) eqs. (5):
             z <= 9.15 m :  r_d = 1.0   - 0.00765 z
             9.15 < z <= 23 m : r_d = 1.174 - 0.0267  z
             23 < z <= 30 m :   r_d = 0.744 - 0.008   z
             z > 30 m :         r_d = 0.5
         (z in meters below ground surface.)

      2) Cyclic Resistance Ratio (CRR_7.5), NCEER SPT clean-sand base curve,
         Youd et al. (2001) eq. (4):
             CRR_7.5 = 1/(34 - (N1)60cs) + (N1)60cs/135
                       + 50/(10*(N1)60cs + 45)^2 - 1/200
         Valid for (N1)60cs < 30; at (N1)60cs >= 30 the soil is treated as too
         dense to liquefy (FS reported as a large sentinel).

      3) Fines correction to (N1)60cs, Youd et al. (2001) eqs. (6a-c, 7a-c):
             (N1)60cs = alpha + beta * (N1)60
         with alpha, beta functions of fines content FC (%).

      4) Magnitude Scaling Factor (MSF), Idriss recommendation in
         Youd et al. (2001) eq. (24):
             MSF = 10^2.24 / M_w^2.56

      5) Overburden correction K_sigma, Hynes-Olsen in Youd et al. (2001):
             K_sigma = (sigma'_v / Pa)^(f-1),  f ~ 0.7-0.8 for typical sands
         Computed only when sigma'_v > Pa; else K_sigma = 1. K_alpha (sloping
         ground) is left at 1.0 (level ground).

      6) Factor of Safety against triggering:
             FS = (CRR_7.5 / CSR) * MSF * K_sigma

    Susceptibility screening uses the soil-type codes in borings_data:
      types 1-4 (sands, gravelly sand) are screened in; type 5 (clayey sand)
      and type 6 (silts) are screened in only if FC < 35% — borderline
      candidates that still need the Bray and Sancio (2006) check the user
      should review. Layers above the groundwater table are unsaturated and
      not analyzed.
    """

    G_FT_S2 = 32.174   # gravity, ft/s^2 (for a_max input in ft/s^2 or as fraction of g)
    PA_PSF = 2116.0    # atmospheric pressure, psf
    M_PER_FT = 0.3048
    FS_NON_LIQUEFIABLE = 99.0   # sentinel for non-susceptible / too-dense layers

    DEFAULT_FINES = {
        1: 5.0,    # Sand (NC)
        2: 5.0,    # Sand (saturated)
        3: 5.0,    # Sand (OC)
        4: 3.0,    # Gravelly Sand
        5: 25.0,   # Clayey Sand
        6: 50.0,   # Silts, sandy silts, clayey silt
    }

    SUSCEPTIBLE_TYPES = {1, 2, 3, 4, 5, 6}

    def __init__(self,
                 borings_data: dict,
                 magnitude: float,
                 peak_acceleration_g: float,
                 gwt_depth: float,
                 unit_weight: float = 120.0,
                 atm_pressure_psf: float = 2116.0,
                 fines_content: Optional[dict] = None,
                 k_sigma_exponent_f: float = 0.75):
        """
        Args:
            borings_data: dict produced by SPT/N60 with keys
                {boring_id: {'depths': [ft], 'n_values': [...],
                             'n60_values': [...], 'n160_values': [...],
                             'soil types': [1..6]}}
            magnitude: design earthquake moment magnitude M_w
            peak_acceleration_g: peak ground surface acceleration as fraction of g
            gwt_depth: groundwater table depth (ft) below ground surface
            unit_weight: representative moist unit weight (pcf) used when a layer
                lacks a unit weight; effective unit weight below the gwt is
                unit_weight - 62.4
            atm_pressure_psf: atmospheric pressure Pa (psf), default 2116
            fines_content: optional dict {boring_id: [FC% per depth]}; if a
                boring or a depth is missing, falls back to DEFAULT_FINES
                keyed by soil type
            k_sigma_exponent_f: exponent f in K_sigma = (sigma'_v/Pa)^(f-1);
                Youd et al. (2001) suggests f ~ 0.7-0.8 for typical sands
        """
        self.borings_data = borings_data
        self.M_w = magnitude
        self.amax_g = peak_acceleration_g
        self.gwt_depth = gwt_depth
        self.gamma = unit_weight
        self.Pa = atm_pressure_psf
        self.fines = fines_content or {}
        self.f_exp = k_sigma_exponent_f

    @staticmethod
    def _rd(z_ft: float) -> float:
        """Liao and Whitman (1986) stress reduction coefficient; z input in ft, converted to m."""
        z_m = z_ft * Liquefaction.M_PER_FT
        if z_m <= 9.15:
            return 1.0 - 0.00765 * z_m
        if z_m <= 23.0:
            return 1.174 - 0.0267 * z_m
        if z_m <= 30.0:
            return 0.744 - 0.008 * z_m
        return 0.5

    def _stresses(self, depth_ft: float) -> Tuple[float, float]:
        """Return (total stress sigma_v, effective stress sigma'_v) at depth, in psf."""
        sigma_v = self.gamma * depth_ft
        if depth_ft <= self.gwt_depth:
            sigma_v_eff = sigma_v
        else:
            u = 62.4 * (depth_ft - self.gwt_depth)
            sigma_v_eff = sigma_v - u
        return sigma_v, sigma_v_eff

    @staticmethod
    def _fines_correction(N1_60: float, FC: float) -> float:
        """Youd et al. (2001) eqs. (6a-c, 7a-c): convert (N1)60 to equivalent clean-sand (N1)60cs."""
        if FC <= 5.0:
            alpha = 0.0
            beta = 1.0
        elif FC < 35.0:
            alpha = math.exp(1.76 - 190.0 / (FC ** 2))
            beta = 0.99 + (FC ** 1.5) / 1000.0
        else:
            alpha = 5.0
            beta = 1.2
        return alpha + beta * N1_60

    @staticmethod
    def _crr_7p5(N1_60cs: float) -> float:
        """NCEER SPT clean-sand base curve, Youd et al. (2001) eq. (4)."""
        if N1_60cs >= 30.0:
            return float('inf')   # treat as non-liquefiable
        return (1.0 / (34.0 - N1_60cs)
                + N1_60cs / 135.0
                + 50.0 / (10.0 * N1_60cs + 45.0) ** 2
                - 1.0 / 200.0)

    def _msf(self) -> float:
        """Idriss MSF per Youd et al. (2001) eq. (24)."""
        return (10.0 ** 2.24) / (self.M_w ** 2.56)

    def _k_sigma(self, sigma_v_eff_psf: float) -> float:
        """Hynes-Olsen K_sigma; capped at 1.0 for sigma'_v <= Pa."""
        if sigma_v_eff_psf <= self.Pa:
            return 1.0
        return (sigma_v_eff_psf / self.Pa) ** (self.f_exp - 1.0)

    def _fc_for(self, boring_id: str, idx: int, soil_type: int) -> float:
        per_boring = self.fines.get(boring_id)
        if per_boring is not None and idx < len(per_boring) and per_boring[idx] is not None:
            return float(per_boring[idx])
        return self.DEFAULT_FINES.get(soil_type, 15.0)

    def analyze_boring(self, boring_id: str) -> dict:
        """Run liquefaction triggering analysis for one boring."""
        if boring_id not in self.borings_data:
            raise ValueError(f"Boring '{boring_id}' not found")
        data = self.borings_data[boring_id]
        depths = data.get('depths')
        n160 = data.get('n160_values')
        soil_types = data.get('soil types') or data.get('soil_types')
        if not depths or not n160 or not soil_types:
            raise ValueError(
                f"Boring '{boring_id}' missing required fields. "
                "Need depths, n160_values, and soil types (run SPT -> N60 first)."
            )

        msf = self._msf()
        rows = []
        for i, (z, N160, st) in enumerate(zip(depths, n160, soil_types)):
            sigma_v, sigma_v_eff = self._stresses(z)
            rd = self._rd(z)

            if z <= self.gwt_depth:
                rows.append({
                    'depth_ft': z, 'soil_type': st, 'N160': N160,
                    'sigma_v_psf': sigma_v, 'sigma_v_eff_psf': sigma_v_eff,
                    'r_d': rd, 'CSR': None, 'CRR_7_5': None,
                    'N1_60cs': None, 'FC_pct': None,
                    'MSF': msf, 'K_sigma': None, 'FS': None,
                    'status': 'above groundwater - not liquefiable',
                })
                continue

            if st not in self.SUSCEPTIBLE_TYPES:
                rows.append({
                    'depth_ft': z, 'soil_type': st, 'N160': N160,
                    'sigma_v_psf': sigma_v, 'sigma_v_eff_psf': sigma_v_eff,
                    'r_d': rd, 'CSR': None, 'CRR_7_5': None,
                    'N1_60cs': None, 'FC_pct': None,
                    'MSF': msf, 'K_sigma': None, 'FS': self.FS_NON_LIQUEFIABLE,
                    'status': f'soil type {st} not susceptible',
                })
                continue

            FC = self._fc_for(boring_id, i, st)

            if st in (5, 6) and FC >= 35.0:
                rows.append({
                    'depth_ft': z, 'soil_type': st, 'N160': N160,
                    'sigma_v_psf': sigma_v, 'sigma_v_eff_psf': sigma_v_eff,
                    'r_d': rd, 'CSR': None, 'CRR_7_5': None,
                    'N1_60cs': None, 'FC_pct': FC,
                    'MSF': msf, 'K_sigma': None, 'FS': self.FS_NON_LIQUEFIABLE,
                    'status': f'fines-dominated (FC={FC:.0f}%) - check Bray & Sancio (2006)',
                })
                continue

            CSR = 0.65 * self.amax_g * (sigma_v / sigma_v_eff) * rd
            N1_60cs = self._fines_correction(N160, FC)
            CRR = self._crr_7p5(N1_60cs)
            K_sig = self._k_sigma(sigma_v_eff)

            if math.isinf(CRR):
                FS = self.FS_NON_LIQUEFIABLE
                status = '(N1)60cs >= 30 - too dense to liquefy'
            elif CSR <= 0:
                FS = self.FS_NON_LIQUEFIABLE
                status = 'CSR <= 0'
            else:
                FS = (CRR / CSR) * msf * K_sig
                status = 'liquefiable' if FS < 1.0 else 'marginal' if FS < 1.3 else 'safe'

            rows.append({
                'depth_ft': z, 'soil_type': st, 'N160': N160,
                'sigma_v_psf': sigma_v, 'sigma_v_eff_psf': sigma_v_eff,
                'r_d': rd, 'CSR': CSR, 'CRR_7_5': CRR if math.isfinite(CRR) else None,
                'N1_60cs': N1_60cs, 'FC_pct': FC,
                'MSF': msf, 'K_sigma': K_sig, 'FS': FS,
                'status': status,
            })

        return {'boring_id': boring_id, 'M_w': self.M_w, 'amax_g': self.amax_g,
                'gwt_depth_ft': self.gwt_depth, 'rows': rows}

    def analyze_all(self) -> dict:
        return {bid: self.analyze_boring(bid) for bid in self.borings_data}

    def run(self):
        """Interactive terminal flow. Expects borings_data populated by SPT and N60."""
        if not self.borings_data:
            print("Error: no borings data. Run SPT -> N60 first.")
            return
        missing = [bid for bid, d in self.borings_data.items() if not d.get('n160_values')]
        if missing:
            print("Error: missing n160_values in borings: " + ", ".join(missing))
            print("Run N60 calculator first.")
            return

        print("\n" + "=" * 70)
        print("Liquefaction Analysis - Youd et al. (2001)")
        print("=" * 70)
        try:
            self.M_w = float(input(f"Earthquake magnitude M_w [{self.M_w}]: ").strip() or self.M_w)
            self.amax_g = float(input(f"Peak ground accel a_max (g) [{self.amax_g}]: ").strip() or self.amax_g)
            self.gwt_depth = float(input(f"Groundwater table depth (ft) [{self.gwt_depth}]: ").strip() or self.gwt_depth)
            self.gamma = float(input(f"Moist unit weight (pcf) [{self.gamma}]: ").strip() or self.gamma)
            self.Pa = float(input(f"Atmospheric pressure (psf) [{self.Pa}]: ").strip() or self.Pa)
            self.f_exp = float(input(f"K_sigma exponent f (0.7-0.8) [{self.f_exp}]: ").strip() or self.f_exp)
        except ValueError as e:
            print(f"Error: invalid input ({e})")
            return

        for bid in self.borings_data:
            ans = input(f"Override fines content (FC %) for {bid}? (y/N): ").strip().lower()
            if ans.startswith('y'):
                n = len(self.borings_data[bid].get('depths', []))
                try:
                    raw = input(f"  Enter {n} FC values (comma-separated, %): ").strip()
                    fc_list = [float(x.strip()) for x in raw.split(',')]
                    if len(fc_list) != n:
                        print(f"  Expected {n} values, got {len(fc_list)}. Using soil-type defaults.")
                    else:
                        self.fines[bid] = fc_list
                except ValueError:
                    print("  Invalid input; using soil-type defaults.")

        self.print_report()

    def print_report(self) -> dict:
        results = self.analyze_all()
        print("\n" + "=" * 110)
        print("LIQUEFACTION ANALYSIS - Youd et al. (2001)")
        print("=" * 110)
        print(f"  Earthquake magnitude M_w : {self.M_w:.2f}")
        print(f"  Peak ground accel a_max  : {self.amax_g:.3f} g")
        print(f"  Groundwater table        : {self.gwt_depth:.2f} ft")
        print(f"  Unit weight (moist)      : {self.gamma:.1f} pcf")
        print(f"  Atmospheric pressure Pa  : {self.Pa:.0f} psf")
        print(f"  K_sigma exponent f       : {self.f_exp:.2f}")

        for bid, r in results.items():
            print("\n" + "-" * 110)
            print(f"BORING: {bid}")
            print("-" * 110)
            header = (f"{'Depth':<7} {'Type':<5} {'N160':<7} "
                      f"{'sigv':<8} {chr(963)+chr(39)+'v':<8} {'rd':<6} "
                      f"{'CSR':<7} {'(N1)cs':<8} {'CRR':<7} {'K_sig':<7} "
                      f"{'FS':<7} {'Status':<32}")
            print(header)
            for row in r['rows']:
                def f(v, fmt):
                    return format(v, fmt) if isinstance(v, (int, float)) and v is not None else 'N/A'
                print(f"{row['depth_ft']:<7.1f} {row['soil_type']:<5} "
                      f"{row['N160']:<7.2f} "
                      f"{row['sigma_v_psf']:<8.0f} {row['sigma_v_eff_psf']:<8.0f} "
                      f"{row['r_d']:<6.3f} "
                      f"{f(row['CSR'], '<7.3f')} "
                      f"{f(row['N1_60cs'], '<8.2f')} "
                      f"{f(row['CRR_7_5'], '<7.3f')} "
                      f"{f(row['K_sigma'], '<7.3f')} "
                      f"{f(row['FS'], '<7.2f')} "
                      f"{row['status']:<32}")
        print("\n" + "=" * 110)
        print("Notes: FS < 1.0 = liquefiable, 1.0-1.3 = marginal, > 1.3 = safe.")
        print("       FS = 99.0 sentinel marks layers screened out (non-susceptible / too dense).")
        print("       Saturated silts/clayey sands with FC >= 35% should also be checked")
        print("       against Bray & Sancio (2006) criteria.")
        print("=" * 110)
        return results


if __name__ == "__main__":
    main_menu()