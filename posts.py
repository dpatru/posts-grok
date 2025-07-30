import math
from math import ceil

def select_hss_size(P_total, L, beam_load, beam_width, F_y=46, phi_c=0.9):
    """
    Select smallest HSS size for given load, length, and eccentricity.
    
    Parameters:
    - P_total: Total factored load (kips)
    - L: Post length (in)
    - beam_load: LVL load causing eccentricity (kips)
    - beam_width: LVL width (in)
    - F_y: HSS yield strength (ksi)
    - phi_c: Compression strength reduction factor
    
    Returns: (HSS designation, area, section modulus, width)
    """
    hss_sizes = [
        ("HSS6x6x3/16", 4.30, 10.3, 6.0, 12.5),
        ("HSS8x8x3/16", 5.75, 17.4, 8.0, 16.8),
        ("HSS8x8x1/4", 7.56, 22.4, 8.0, 22.0),
        ("HSS10x10x1/4", 9.52, 35.2, 10.0, 27.5),
        ("HSS12x12x1/4", 11.5, 50.8, 12.0, 33.0)
    ]
    
    E = 29000  # Steel modulus (ksi)
    K = 1.0    # Effective length factor
    r_factor = 0.4  # Approx radius of gyration factor
    
    for hss, A, S, B_hss, _ in hss_sizes:
        phi_Pn = phi_c * F_y * A
        slenderness = K * L / (r_factor * math.sqrt(A))
        if slenderness > 200:
            continue
        # Eccentricity: LVL edge at HSS edge
        e = B_hss / 2 - beam_width / 2
        M_u = beam_load * e
        phi_Mn = phi_c * F_y * S
        interaction = (P_total / phi_Pn) + (M_u / phi_Mn)
        # Add safety margin: interaction <= 0.9
        if interaction <= 0.9:
            return hss, A, S, B_hss
    return None, None, None, None

def calculate_plates(P_beam, P_upper, B_hss, beam_width, beam_height, beam_offset=0, 
                    lvl_bearing=500, F_y=36, phi=0.9, load_factor=1.6):
    """
    Calculate cap and bearing plates for one stack level.
    
    Parameters:
    - P_beam: Factored LVL load (kips)
    - P_upper: Factored load from above (kips)
    - B_hss: HSS width (in)
    - beam_width: LVL width (in)
    - beam_height: LVL depth (in)
    - beam_offset: Distance from HSS edge to LVL edge (in)
    - lvl_bearing: LVL bearing stress (psi)
    - F_y: Plate yield strength (ksi)
    - phi: Strength reduction factor
    - load_factor: LRFD factor
    
    Returns: dict with plate dimensions and stresses
    """
    try:
        # LVL bearing area
        P_s = P_beam / load_factor
        A_bearing = P_s * 1000 / lvl_bearing
        N = A_bearing / beam_width
        N = max(N, B_hss)
        N = ceil(N * 4) / 4  # Round up to ensure bearing stress <= 500 psi
        
        # Cap plate width
        B = max(B_hss + 2, beam_width + 2)
        
        # Cap plate thickness
        d = B_hss
        m = (B - 0.95 * d) / 2
        n = (B - 0.8 * beam_width) / 2
        beam_centroid = beam_offset + beam_width / 2
        n_ecc = max(B_hss - (beam_offset + beam_width), beam_centroid)
        l = max(m, n, n_ecc)
        
        denominator = phi * F_y * B * N
        if denominator == 0:
            raise ValueError("Invalid plate dimensions.")
        t_cap = l * ((2 * (P_beam + P_upper)) / denominator) ** 0.5
        t_cap = max(round(t_cap * 8) / 8, 0.375)  # Minimum 3/8" thickness
        if t_cap < l * ((2 * (P_beam + P_upper)) / denominator) ** 0.5:
            t_cap += 0.125
        
        # Bearing plates
        bearing_width = B_hss
        bearing_height = beam_height
        t_bearing = 1.0
        bearing_area = bearing_width * t_bearing
        phi_Rn = 0.75 * 1.8 * F_y * bearing_area
        if phi_Rn < P_upper / 2:
            t_bearing = (P_upper / 2) / (0.75 * 1.8 * F_y * bearing_width)
            t_bearing = round(t_bearing * 8) / 8
            if t_bearing < (P_upper / 2) / (0.75 * 1.8 * F_y * bearing_width):
                t_bearing += 0.125
        
        # LVL bearing stress
        bearing_stress = P_s * 1000 / (beam_width * N)
        if bearing_stress > lvl_bearing:
            raise ValueError(f"Bearing stress {bearing_stress:.1f} psi exceeds {lvl_bearing} psi.")
        
        # Welds
        weld_size = 3/16
        phi_Rn_weld = 0.75 * 0.707 * weld_size * 0.6 * 70
        q_cap = (P_beam + P_upper) / (4 * B_hss)
        q_bearing = P_upper / (2 * 2 * B_hss)
        if max(q_cap, q_bearing) > phi_Rn_weld:
            weld_size = max(q_cap, q_bearing) / (0.75 * 0.707 * 0.6 * 70)
            weld_size = round(weld_size * 16) / 16
        
        return {
            "Cap Plate": f"{B:.3f} x {N:.3f} x {t_cap:.3f} in",
            "Bearing Plates": f"2 plates, {bearing_width:.3f} x {bearing_height:.3f} x {t_bearing:.3f} in",
            "LVL Bearing Stress": f"{bearing_stress:.1f} psi",
            "Weld Size": f"{weld_size:.3f} in"
        }
    
    except Exception as e:
        return {"Error": str(e)}

def process_column_stack(cases, lvl_bearing=500):
    """
    Process column-beam stack for HSS, cap plates, and bearing plates.
    
    Parameters:
    - cases: List of (post_length, beam_load, beam_width, beam_height)
    - lvl_bearing: LVL bearing stress (psi)
    
    Returns: List of design results
    """
    results = []
    # Calculate cumulative loads
    total_loads = []
    cumulative_load = 0
    for _, beam_load, _, _ in reversed(cases):
        cumulative_load += beam_load
        total_loads.insert(0, cumulative_load)
    
    # Process each stack level
    for i, ((post_length, beam_load, beam_width, beam_height), P_total) in enumerate(zip(cases, total_loads), 1):
        # Select HSS
        hss, A_hss, S_hss, B_hss = select_hss_size(P_total, post_length, beam_load, beam_width)
        if hss is None:
            results.append(f"Level {i}: {{'Error': 'No suitable HSS size found'}}")
            continue
        
        # Calculate plates (P_upper = load from all levels above)
        P_upper = P_total - beam_load
        plate_design = calculate_plates(beam_load, P_upper, B_hss, beam_width, beam_height, 
                                       beam_offset=0, lvl_bearing=lvl_bearing)
        
        results.append(f"Level {i}: {{'HSS': '{hss}', {', '.join(f'\'{k}\': \'{v}\'' for k, v in plate_design.items())}}}")
    
    return results

# Example usage
if __name__ == "__main__":
    # Input: (post_length, beam_load, beam_width, beam_height)
    cases = [
        (120, 56, 3, 12),    # Level 1
        (96, 40, 3.5, 14),   # Level 2
        (144, 20, 2.5, 10),  # Level 3
        (144, 0.1, 2.5, 10)  # Level 4
    ]
    
    results = process_column_stack(cases)
    for result in results:
        print(result)
        print()
